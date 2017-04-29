#pragma once

#include <itkImageToImageFilter.h>

using namespace agtk;

template<typename TInputImage, typename TOutputImage>
class AirwaysMeasureImageFilter : public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  typedef AirwaysMeasureImageFilter Self;

  typedef itk::ImageToImageFilter<TInputImage, TOutputImage> Superclass;

  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;

  typedef typename InputImageType::PixelType   InputPixelType;
  typedef typename OutputImageType::PixelType  OutputPixelType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef double EigenValueType;
  typedef itk::FixedArray<EigenValueType, InputImageType::ImageDimension> EigenValueArrayType;
  typedef itk::SymmetricEigenAnalysis<InputPixelType, EigenValueArrayType> EigenCalculatorType;

  typedef itk::ImageRegionConstIterator<InputImageType> InputImageRegionIterator;
  typedef itk::ImageRegionIterator<OutputImageType> OutputImageRegionIterator;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(AirwaysMeasureImageFilter, ImageToImageFilter);

  /** Image dimension */
  itkStaticConstMacro(ImageDimension, unsigned int, InputImageType::ImageDimension);

  static_assert(ImageDimension == 3U, "Invalid image dimension. 3D images are supported.");

protected:
  AirwaysMeasureImageFilter();
  ~AirwaysMeasureImageFilter() {}

  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
  AirwaysMeasureImageFilter(const Self &);  // purposely not implemented
  void operator=(const Self &);                    // purposely not implemented

  double ComputeTubeValue(const InputPixelType& pixelValue);
  void ComputeEigenValues(const InputPixelType& A, EigenValueArrayType& eigenValues);

  double ComputePhiValue(double e2, double e3);
  double ComputeWeightValue(double e1, double e2);

  double m_Mu = 1;
  double m_Alpha = 0.25;
};
#ifndef ITK_MANUAL_INSTANTIATION
#include "AirwaysMeasureFilter.hxx"
#endif