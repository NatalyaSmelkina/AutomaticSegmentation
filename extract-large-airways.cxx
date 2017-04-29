#include <itkLaplacianSharpeningImageFilter.h>
#include <itkTimeProbe.h>
#include <itkLaplacianRecursiveGaussianImageFilter.h>
#include <itkImageMomentsCalculator.h>
#include <itkConnectedThresholdImageFilter.h>

#include <agtkTypes.h>
#include <agtkIO.h>
#include <agtkPath.h>
#include <agtkCommandLineArgumentParser.h>
#include <agtkBinaryImageUtilities.h>

#include "AirwaysMeasureFilter.h"
#include "HessianAnalysisFilter.h"

using namespace agtk;

int main(int argc, char* argv[])
{
  CommandLineArgumentParser::Pointer parser = agtk::CommandLineArgumentParser::New();
  parser->SetCommandLineArguments(argc, argv);

  std::string inputFilename;
  parser->GetValue("-image", inputFilename);

  std::string outputFilename;
  parser->GetValue("-output", outputFilename);

  /*read image*/
  FloatImage3D::Pointer image = FloatImage3D::New();
  if (!readImage<FloatImage3D>(image, inputFilename)) {
    return EXIT_FAILURE;
  }

  itk::TimeProbe clock;
  clock.Start();
  /*sharp image*/
  /*typedef itk::LaplacianSharpeningImageFilter<FloatImage3D, FloatImage3D >  LaplacianSharpeningImageFilterType;
  LaplacianSharpeningImageFilterType::Pointer laplacianSharpeningImageFilter =
    LaplacianSharpeningImageFilterType::New();

  laplacianSharpeningImageFilter->SetInput(image);
  */
  /*if (!writeImage<FloatImage3D>(laplacianSharpeningImageFilter->GetOutput(), addFileNameSuffix(outputFilename, "-sharpened"))) {
    return EXIT_FAILURE;
  }*/

  FloatImage3D::IndexType index;
  index[0] = 0;
  index[1] = 0;
  index[2] = 0;
  auto weight = image->GetLargestPossibleRegion().GetSize()[0];
  auto hight = image->GetLargestPossibleRegion().GetSize()[1];
  auto length = image->GetLargestPossibleRegion().GetSize()[2];

  typedef itk::LaplacianRecursiveGaussianImageFilter<FloatImage3D, FloatImage3D>  LaplacianRecursiveGaussianImageFilterType;
  LaplacianRecursiveGaussianImageFilterType::Pointer LoGFilter = LaplacianRecursiveGaussianImageFilterType::New();
  LoGFilter->SetInput(image);
  LoGFilter->SetSigma(0.5/image->GetSpacing()[0]);
  LoGFilter->SetNormalizeAcrossScale(false);
  LoGFilter->Update();
  FloatImage3D::Pointer LoGImage = LoGFilter->GetOutput();
  float w = 0.0;
  for (auto k = 0; k < length; k++) {
    for (auto j = 0; j < hight; j++) {
      for (auto i = 0; i < weight; i++) {
        index[0] = i;
        index[1] = j;
        index[2] = k;
        if (LoGImage->GetPixel(index) <= 0) {
          w = 0.5;
        }
        else {
          w = 0.05;
        }
        LoGImage->SetPixel(index, image->GetPixel(index) - w*LoGImage->GetPixel(index));
      }
    }
  }
  /*if (!writeImage<FloatImage3D>(LoGImage, addFileNameSuffix(outputFilename, "-sharpened-by-LoG"))) {
    return EXIT_FAILURE;
  }*/

  HessianAnalysisImageFilter::Pointer airwayEnhancer = HessianAnalysisImageFilter::New();
  airwayEnhancer->SetInput(LoGImage);

  try {
    airwayEnhancer->Update();
  }
  catch (itk::ExceptionObject&) {
    return EXIT_FAILURE;
  }

  BinaryImage3D::Pointer trachea = binarizeImage(airwayEnhancer->GetOutput(), 125, FloatLimits::infinity()); // to do low treshold
  trachea = getLargestObjectFromBinaryImage(trachea);

  if (!writeImage<BinaryImage3D>(trachea, addFileNameSuffix(outputFilename, "-enhanced-trachea-by-log"))) {
    return EXIT_FAILURE;
  }

  typedef itk::ImageMomentsCalculator< BinaryImage3D > ImageCalculatorType;
  ImageCalculatorType::Pointer airwayCalculator = ImageCalculatorType::New();
  airwayCalculator->SetImage(trachea);
  airwayCalculator->Compute();
  ImageCalculatorType::VectorType center = airwayCalculator->GetCenterOfGravity();

  BinaryImage3D::PointType point;
  point[0] = center[0];
  point[1] = center[1];
  point[2] = center[2];
  BinaryImage3D::IndexType pixelIndex;
  bool transform = trachea->TransformPhysicalPointToIndex(point, pixelIndex);

  typedef itk::ConnectedThresholdImageFilter< FloatImage3D, BinaryImage3D > ConnectedFilterType;
  ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
  connectedThreshold->SetInput(image);
  connectedThreshold->SetLower(-FloatLimits::infinity());
  connectedThreshold->SetUpper(-600);
  connectedThreshold->SetSeed(pixelIndex);
  connectedThreshold->SetReplaceValue(255);

  if (!writeImage<BinaryImage3D>(connectedThreshold->GetOutput(), addFileNameSuffix(outputFilename, "-init-lungs-by-log"))) {
    return EXIT_FAILURE;
  }

  clock.Stop();
  std::cout << "elapsed time, sec " << clock.GetTotal() << std::endl;

  return EXIT_SUCCESS;
}