#pragma once

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkSymmetricSecondRankTensor.h>
#include <itkMultiScaleHessianBasedMeasureImageFilter.h>

#include "AirwaysMeasureFilter.h"
#include "HessianAnalysisFilter.h"

HessianAnalysisImageFilter::HessianAnalysisImageFilter()
  : m_SigmaMinimum(0.5)
  , m_SigmaMaximum(5)
  , m_SigmaSteps(10)
  , m_ProgressAccumulator(itk::ProgressAccumulator::New())
{
}

void HessianAnalysisImageFilter::GenerateData()
{
  m_ProgressAccumulator->ResetProgress();
  m_ProgressAccumulator->SetMiniPipelineFilter(this);
  m_ProgressAccumulator->UnregisterAllFilters();

  ComputeTubeLikeStructure();
}

void HessianAnalysisImageFilter::ComputeTubeLikeStructure()
{
  typedef itk::SymmetricSecondRankTensor<double, IMAGE_DIM_3> HessianPixel3D;
  typedef itk::Image<HessianPixel3D, IMAGE_DIM_3> HessianImage3D;

  typedef AirwaysMeasureImageFilter<HessianImage3D, FloatImage3D> AirwaysMeasurer;
  auto airwaysMeasurer = AirwaysMeasurer::New();

  typedef itk::MultiScaleHessianBasedMeasureImageFilter<FloatImage3D, HessianImage3D, FloatImage3D> MultiScaleHessianBasedMeasureImageFilter;
  auto multiScaleMeasurer = MultiScaleHessianBasedMeasureImageFilter::New();

  m_ProgressAccumulator->RegisterInternalFilter(multiScaleMeasurer, 1.0);

  multiScaleMeasurer->SetHessianToMeasureFilter(airwaysMeasurer);
  multiScaleMeasurer->SetSigmaMinimum(m_SigmaMinimum);
  multiScaleMeasurer->SetSigmaMaximum(m_SigmaMaximum);
  multiScaleMeasurer->SetNumberOfSigmaSteps(m_SigmaSteps);

  multiScaleMeasurer->SetSigmaStepMethodToEquispaced();

  multiScaleMeasurer->SetInput(this->GetInput());
  multiScaleMeasurer->Update();

  this->GraftOutput(multiScaleMeasurer->GetOutput());
}