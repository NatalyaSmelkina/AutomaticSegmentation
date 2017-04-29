#pragma once

#include <itkProgressAccumulator.h>

#include <itkImageToImageFilter.h>

using namespace agtk;

class HessianAnalysisImageFilter : public itk::ImageToImageFilter<FloatImage3D, FloatImage3D>
{
public:
  typedef HessianAnalysisImageFilter Self;
  typedef itk::ImageToImageFilter<FloatImage3D, FloatImage3D> Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  itkNewMacro(Self);
  itkTypeMacro(HessianAnalysisImageFilter, ImageToImageFilter);

  itkSetMacro(SigmaMinimum, double);
  itkGetConstMacro(SigmaMinimum, double);

  itkSetMacro(SigmaMaximum, double);
  itkGetConstMacro(SigmaMaximum, double);

  itkSetMacro(SigmaSteps, unsigned int);
  itkGetConstMacro(SigmaSteps, unsigned int);

protected:
  HessianAnalysisImageFilter();
  virtual ~HessianAnalysisImageFilter() {}

  virtual void GenerateData() override;

  double m_SigmaMinimum;
  double m_SigmaMaximum;
  unsigned int m_SigmaSteps;
  itk::ProgressAccumulator::Pointer m_ProgressAccumulator;

private:
  HessianAnalysisImageFilter(const Self&);
  void operator=(const Self&);

  void ComputeTubeLikeStructure();
};
#ifndef ITK_MANUAL_INSTANTIATION
#include "HessianAnalysisFilter.hxx"
#endif