SET (src_dir ${CMAKE_CURRENT_SOURCE_DIR})

SET (tpc_phase_shift_calibrator_sources
  ${src_dir}/TPCPhaseShiftCalibratorAL.cc
  ${src_dir}/TPCPhaseShiftCalibratorAL.h
  PARENT_SCOPE
)

SET (tpc_phase_shift_calibrator_data
  ${src_dir}/TPCPhaseShiftCalibratorAL.xml
  ${src_dir}/TPCPhaseShiftCalibratorAL.xsd
  PARENT_SCOPE
)
