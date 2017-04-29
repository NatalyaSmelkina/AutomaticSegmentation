#pragma once

#include <cmath>

#include <itkProgressReporter.h>

#include "AirwaysMeasureFilter.h"

#define SQRT3 1.73205080756887729352744634151  // sqrt(3)
#define SQR(x) ((x)*(x))

template<typename TInputImage, typename TOutputImage>
AirwaysMeasureImageFilter<TInputImage, TOutputImage>
::AirwaysMeasureImageFilter()
{
}

template<typename TInputImage, typename TOutputImage>
double AirwaysMeasureImageFilter<TInputImage, TOutputImage>
::ComputeTubeValue(const InputPixelType& pixelValue)
{
  EigenValueArrayType eigenValues;
  ComputeEigenValues(pixelValue, eigenValues);

  EigenValueType e1 = eigenValues[0];
  EigenValueType e2 = eigenValues[1];
  EigenValueType e3 = eigenValues[2];
  double tubeValue = 0;
  if ((e3 >= e2) && (e2 > 0)) {
    tubeValue = std::abs(e3)*ComputePhiValue(e2, e3);// *ComputeWeightValue(e1, e2);
  }
  return tubeValue;
}

template<typename TInputImage, typename TOutputImage>
double AirwaysMeasureImageFilter<TInputImage, TOutputImage>
::ComputePhiValue(double e2, double e3)
{
  double phi = 0;
  if ((e3 >= e2) && (e2 > 0)) {
    phi = std::pow(e2 / e3, m_Mu);
  }
  return phi;
}

template<typename TInputImage, typename TOutputImage>
double AirwaysMeasureImageFilter<TInputImage, TOutputImage>
::ComputeWeightValue(double e1, double e2)
{
  double w = 0;
  double e = e1 / std::abs(e2);
  if ((e2 >= e1) && (e1 > 0)) {
    w = std::pow((1 + e), m_Mu);
  }
  else if ((e1 < 0) && (e1 > std::abs(e2) / m_Alpha)) {
    w = std::pow((1 - m_Alpha*e), m_Mu);
  }
  return w;
}

template<typename TInputImage, typename TOutputImage>
void AirwaysMeasureImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, itk::ThreadIdType threadId)
{
  typename InputImageType::ConstPointer input = this->GetInput();
  typename OutputImageType::Pointer output = this->GetOutput();

  itk::ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels(), 1000 / this->GetNumberOfThreads());

  InputImageRegionIterator iit(input, outputRegionForThread);
  OutputImageRegionIterator oit(output, outputRegionForThread);

  iit.GoToBegin();
  oit.GoToBegin();

    while (!iit.IsAtEnd()) {
      double tubeValue = ComputeTubeValue(iit.Get());
      oit.Set(tubeValue);

      ++iit;
      ++oit;

      progress.CompletedPixel();
    }
}

template<typename TInputImage, typename TOutputImage>
void AirwaysMeasureImageFilter<TInputImage, TOutputImage>
::ComputeEigenValues(const InputPixelType& A, EigenValueArrayType& w)
{
  //EigenCalculatorType symmetricEigenSystem(InputImageType::ImageDimension);
  //symmetricEigenSystem.ComputeEigenValuesAndVectors(A, w, v);
  double m, c1, c0;

  // Determine coefficients of characteristic poynomial. We write
  //       | a   d   f  |
  //  A =  | d*  b   e  |
  //       | f*  e*  c  |
  double de = A(0, 1) * A(1, 2);  // d * e
  double dd = SQR(A(0, 1));       // d^2
  double ee = SQR(A(1, 2));       // e^2
  double ff = SQR(A(0, 2));       // f^2

  m = A(0, 0) + A(1, 1) + A(2, 2);

  // a*b + a*c + b*c - d^2 - e^2 - f^2
  c1 = (A(0, 0) * A(1, 1) + A(0, 0) * A(2, 2) + A(1, 1) * A(2, 2)) - (dd + ee + ff);

  // c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)
  c0 = A(2, 2) * dd + A(0, 0) * ee + A(1, 1) * ff - A(0, 0) * A(1, 1) * A(2, 2) - 2.0 * A(0, 2) * de;

  double p, sqrt_p, q, c, s, phi;

  p = SQR(m) - 3.0*c1;
  q = m*(p - (3.0 / 2.0)*c1) - (27.0 / 2.0)*c0;
  sqrt_p = std::sqrt(std::fabs(p));

  phi = 27.0 * (0.25*SQR(c1)*(p - c1) + c0*(q + 27.0 / 4.0*c0));
  phi = (1.0 / 3.0) * std::atan2(std::sqrt(std::fabs(phi)), q);

  c = sqrt_p*std::cos(phi);
  s = (1.0 / SQRT3)*sqrt_p*std::sin(phi);

  w[1] = (1.0 / 3.0)*(m - c);
  w[2] = w[1] + s;
  w[0] = w[1] + c;
  w[1] -= s;

  std::sort(w.Begin(), w.End());
}