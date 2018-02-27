SET (src_dir ${CMAKE_CURRENT_SOURCE_DIR})

SET (mhtdc_calibrator_sources
  ${src_dir}/MHTDCCalibratorAL.cc
  ${src_dir}/MHTDCCalibratorAL.h
  PARENT_SCOPE
)

SET (mhtdc_calibrator_data
  ${src_dir}/MHTDCCalibratorAL.xml
  ${src_dir}/MHTDCCalibratorAL.xsd
  PARENT_SCOPE
)
