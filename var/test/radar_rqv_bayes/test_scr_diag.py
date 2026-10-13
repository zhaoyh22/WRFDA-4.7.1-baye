"""Compile the production SCR writer and inspect real NetCDF output."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[2]
RADAR = ROOT / "da/da_radar"
SOURCE = RADAR / "da_write_radar_rqv_scr_diag.inc"


class ScrDiagnostics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.work = Path(cls.tmp.name)
        writer = SOURCE.read_text().split(
            "subroutine da_write_radar_rqv_scr_diag_nc(", 1)[1]
        writer = "subroutine da_write_radar_rqv_scr_diag_nc(" + writer
        host = """
module writer_host
implicit none
integer, parameter :: ids=1, ide=2, jds=1, jde=2, kds=1, kde=2
integer, parameter :: filename_len=256, stdout=6
integer :: num_fgat_time
character(len=19) :: analysis_date='2026-10-13_12:00:00'
contains
""" + writer + """
end module
program regression
use writer_host
implicit none
integer :: slot, count
character(len=32) :: arg
real :: valid(2,2,2), candidate(2,2,2), assimilated(2,2,2)
call get_command_argument(1,arg)
read(arg,*) num_fgat_time
! Three representative slots suffice even for a very large slot count.
count=min(num_fgat_time,3)
do slot=count,1,-1
   valid=0.0
   candidate=0.0
   assimilated=0.0
   valid(slot,1,1)=1.0
   candidate(1,slot,1)=1.0
   assimilated(1,1,slot)=1.0
   call da_write_radar_rqv_scr_diag_nc(slot,valid,candidate,assimilated)
end do
end program
"""
        # Use dimension 3 for distinct slot-specific cells in each mask.
        host = host.replace("ide=2", "ide=3").replace("jde=2", "jde=3")
        host = host.replace("kde=2", "kde=3").replace("(2,2,2)", "(3,3,3)")
        (cls.work / "writer.f90").write_text(host)
        flags = shlex.split(subprocess.check_output(
            ["nf-config", "--fflags"], text=True))
        libs = shlex.split(subprocess.check_output(
            ["nf-config", "--flibs"], text=True))
        subprocess.run([os.environ.get("FC", "gfortran"), "-fcheck=all",
                        *flags, str(cls.work / "writer.f90"), *libs,
                        "-o", str(cls.work / "writer")], check=True,
                       cwd=cls.work)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_slots(self, total):
        work = self.work / str(total)
        work.mkdir()
        subprocess.run([str(self.work / "writer"), str(total)],
                       cwd=work, check=True)
        return work

    def check_file(self, path, slot):
        with netcdf_file(path, "r", mmap=False) as nc:
            self.assertEqual(int(nc.time_index.decode().strip()), slot)
            self.assertEqual(nc.analysis_date, b"2026-10-13_12:00:00")
            self.assertEqual(set(nc.variables),
                             {"scr_valid", "scr_candidate", "scr_assimilated"})
            # NetCDF Fortran dimensions appear reversed in the Python reader.
            for name, axis in [("scr_valid", 2), ("scr_candidate", 1),
                               ("scr_assimilated", 0)]:
                expected = np.zeros((3, 3, 3), dtype=np.int8)
                index = [0, 0, 0]
                index[axis] = slot - 1
                expected[tuple(index)] = 1
                self.assertEqual(nc.variables[name].typecode(), "b")
                np.testing.assert_array_equal(nc.variables[name][:], expected)

    def test_single_slot_legacy_name(self):
        work = self.run_slots(1)
        self.assertEqual({p.name for p in work.iterdir()}, {"scr_diag.nc"})
        self.check_file(work / "scr_diag.nc", 1)

    def test_reverse_order_slots_survive(self):
        work = self.run_slots(3)
        self.assertEqual({p.name for p in work.iterdir()},
                         {f"scr_diag_t{i}.nc" for i in (1, 2, 3)})
        for slot in (1, 2, 3):
            self.check_file(work / f"scr_diag_t{slot}.nc", slot)

    def test_large_index_no_format_overflow(self):
        # Exercise the production writer directly with a large slot index.
        source = (self.work / "writer.f90").read_text()
        source = source.replace("count=min(num_fgat_time,3)", "count=1")
        source = source.replace("call da_write_radar_rqv_scr_diag_nc(slot,",
                                "call da_write_radar_rqv_scr_diag_nc(num_fgat_time,")
        test_source = self.work / "large.f90"
        test_source.write_text(source)
        flags = shlex.split(subprocess.check_output(["nf-config", "--fflags"], text=True))
        libs = shlex.split(subprocess.check_output(["nf-config", "--flibs"], text=True))
        subprocess.run([os.environ.get("FC", "gfortran"), "-fcheck=all",
                        *flags, str(test_source), *libs, "-o", str(self.work / "large")],
                       check=True, cwd=self.work)
        work = self.work / "large_output"
        work.mkdir()
        subprocess.run([str(self.work / "large"), "100000"],
                       check=True, cwd=work)
        self.assertEqual({p.name for p in work.iterdir()}, {"scr_diag_t100000.nc"})
        with netcdf_file(work / "scr_diag_t100000.nc", "r", mmap=False) as nc:
            self.assertEqual(int(nc.time_index.decode().strip()), 100000)

    def test_slot_wiring(self):
        driver = (ROOT / "da/da_minimisation/da_get_innov_vector.inc").read_text()
        self.assertIn("do n= num_fgat_time , 1, -1", driver)
        self.assertLess(driver.index("iv%time = n"),
                        driver.index("call da_get_innov_vector_radar"))
        self.assertIn("plocal(iv%time-1)", driver)
        self.assertIn("plocal(iv%time)", driver)
        self.assertIn("da_write_radar_rqv_scr_diag_nc(iv%time", SOURCE.read_text())
        self.assertIn("num_fgat_time", (RADAR / "da_radar.f90").read_text().split("contains")[0])


if __name__ == "__main__":
    unittest.main()
