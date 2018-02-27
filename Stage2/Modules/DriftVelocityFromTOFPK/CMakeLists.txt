SET (src_dir ${CMAKE_CURRENT_SOURCE_DIR})

SET (tpc_drift_velocity_calibrator_sources
  ${src_dir}/DriftVelocityFromTOFPK.cc
  ${src_dir}/DriftVelocityFromTOFPK.h
  PARENT_SCOPE
)

SET (tpc_drift_velocity_calibrator_data
  ${src_dir}/DriftVelocityFromTOFPK.xml
  ${src_dir}/DriftVelocityFromTOFPK.xsd
  PARENT_SCOPE
)
