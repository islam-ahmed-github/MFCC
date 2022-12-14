
 
 
 
 /*
 * MFCC.h
 *
 *  Created on: Mar 20, 2021
 *      Author: Islam Ahmed
 */

 
 
#include "MFCC.h"





MFCC::MFCC()
{
	//ctor
}

MFCC::~MFCC()
{
	//dtor
}



 
 

double GetMelFromFhz(double Fhz)
{
	return  (1125.0 * (std::log(1.0 + (Fhz / 700.0))));
}


double GetFhzFromMel(double Mel)
{
	return (700.0 * (std::exp((Mel) / 1125.0) - 1.0));
}


double GetVectorAbsoluteSum(double *Vector, UINT32 Len)
{
	double sum = 0.0;
	UINT32 i = 0;
	for (i = 0; i<Len; i++)  sum += std::abs(Vector[i]);
	return sum;
}



void MFCC::GetDCT(const double *InputFrame,  UINT32 InputFrameLength, double *OutputFrame  ,  UINT32 NumOfCoeff)
{
	double sum = 0.0;
	UINT32 Cindex, SampleIndex;
	if ((NumOfCoeff <= InputFrameLength) && (InputFrameLength >= 1))
	{
		//Coefficients index
		for (Cindex = 0; Cindex < NumOfCoeff; Cindex++)// iterate over coefficients to be calculated
		{
			sum = 0.0;
			//(Filter Energies Sum) index = sample of the input frame
			for (SampleIndex = 0; SampleIndex < InputFrameLength; SampleIndex++)//iterate over frame samples
			{
				sum += (InputFrame[SampleIndex] * std::cos(((double)Cindex*((double)SampleIndex - 0.5)*MACRO_PI) / (InputFrameLength)));
			}//end of inner for
			OutputFrame[Cindex] = sum;
		}//end outter for
	}
}// end function "GetDCT"



// filter bank width---> 
void MFCC::GetMFCC(const double *InputFrame, const uint32_t InputFrameLength, const double *FiltersBank, UINT32 NumOfFilters , double *OutputFrame)
{
	if ((InputFrameLength>1))
	{
		uint32_t FilterIndex, SampleIndex;
		double sum = 0.0;
		for (FilterIndex = 0; FilterIndex < NumOfFilters; FilterIndex++)// iterate samples of each filters
		{
			sum = 0.0;
			for (SampleIndex = 0; SampleIndex < InputFrameLength; SampleIndex++)// iterate samples of each filters
			{
				sum += (std::abs(InputFrame[SampleIndex] * FiltersBank[MACRO_GetArrayIndex(FilterIndex, SampleIndex, InputFrameLength)]));
				
			}//end for
			OutputFrame[FilterIndex] = std::log10(sum);
		} // end otter for
	}// end if
}// end func




void MFCC::GetMFCCDCT(const double *InputFrame, const uint32_t InputFrameLength, const double *FiltersBank, UINT32 NumOfFilters, UINT32 NumOfCoefficients , double *OutputFrame)
{

	double * temp = new double[NumOfFilters];

	GetMFCC(InputFrame, InputFrameLength, FiltersBank, NumOfFilters, temp);
	GetDCT(temp, NumOfFilters, OutputFrame, NumOfCoefficients);

	delete temp;
}

 


void MFCC::GetDeltaDelta(const double *InputMfcc, uint32_t NumOfVectors, uint32_t VectorLength, double *OutputMfccDD )
{
	 
	if ((NumOfVectors>0) && (VectorLength>0))
	{
		uint32_t  frameIndex, DeltaIndex;
		uint32_t DDOffset = (VectorLength * 2);

		// copy mfcc to the new buffer will apppend d and DD
		for (frameIndex = 0; frameIndex< NumOfVectors; frameIndex++)
		{
			for (DeltaIndex = 0; DeltaIndex <VectorLength; DeltaIndex++)
			{
				OutputMfccDD[MACRO_GetArrayIndex(frameIndex, DeltaIndex, VectorLength*3)] = \
					InputMfcc[MACRO_GetArrayIndex(frameIndex, DeltaIndex, VectorLength)];
			}
		}

		//f0
		for (DeltaIndex = 0; DeltaIndex < VectorLength; DeltaIndex++)
		{
			// delta
			OutputMfccDD[MACRO_GetArrayIndex(0, (VectorLength + DeltaIndex), VectorLength * 3)] = \
				OutputMfccDD[MACRO_GetArrayIndex(0,DeltaIndex, VectorLength * 3)]/2;
			// delta-delta
			OutputMfccDD[MACRO_GetArrayIndex(0, (DDOffset + DeltaIndex), VectorLength * 3)]  = \
				OutputMfccDD[MACRO_GetArrayIndex(0,(VectorLength+DeltaIndex),VectorLength * 3)] / 2;
		}
		//f1--->fn
		for (frameIndex = 1; frameIndex < NumOfVectors; frameIndex++)
		{
			for (DeltaIndex = 0; DeltaIndex < VectorLength; DeltaIndex++)
			{
				// Delta
				OutputMfccDD[MACRO_GetArrayIndex(frameIndex, (VectorLength + DeltaIndex), VectorLength * 3)] = \
					((OutputMfccDD[MACRO_GetArrayIndex(frameIndex, DeltaIndex, VectorLength * 3)] -\
						OutputMfccDD[MACRO_GetArrayIndex(frameIndex - 1, DeltaIndex, VectorLength * 3)]) / 2);

				// Delta Delta
			 	OutputMfccDD[MACRO_GetArrayIndex(frameIndex, (DDOffset + DeltaIndex), VectorLength * 3)] =\
					((OutputMfccDD[MACRO_GetArrayIndex(frameIndex, (VectorLength + DeltaIndex), VectorLength * 3)] -\
						OutputMfccDD[MACRO_GetArrayIndex(frameIndex - 1, (VectorLength + DeltaIndex), VectorLength * 3)]) / 2);
			
			} // end inner for
		} // end outter for
	}// end if
}// end func



 

void MFCC::GenerateFiltersBank(double *Output , UINT32 NumOfFilters , UINT32 SamplingRateHz, UINT32 UpperFreqHz , UINT32 LowerFreqHz , UINT32 FFT_FullSize)
{
	// check arguments validity
	if ((LowerFreqHz>1) && (UpperFreqHz>LowerFreqHz) && (UpperFreqHz <= (SamplingRateHz / 2)) && (FFT_FullSize>2) && (NumOfFilters>2))
	{
		UINT32 FilterIndex, i , FilterLength = ((FFT_FullSize / 2) + 1);

	 
		double  *FilterBankMel = new double[( NumOfFilters + 2)]{ 0 };      // filter bank in  mel scale
		double  *FilterBankHz = new double[( NumOfFilters + 2)]{ 0 };    // filter bank in  frequency scale
		double  *FilterBankFFTbins = new double[(NumOfFilters + 2)]{ 0 };    // filter bank in  FFT-Bins scale


		FilterBankHz[0] = (double)LowerFreqHz; // compute lower Hz
		FilterBankMel[0] = GetMelFromFhz((double)LowerFreqHz);//compute lower mel

		FilterBankHz[ NumOfFilters + 1] = (double)UpperFreqHz;//compute upper Hz
		FilterBankMel[NumOfFilters + 1] = GetMelFromFhz((double)UpperFreqHz);// compute upper mel


		double MelSpace = ((FilterBankMel[ NumOfFilters + 1] - FilterBankMel[0]) / ( NumOfFilters + 1));

		//compute other items in mel and frequency scale
		for (i = 1; i < (NumOfFilters + 1); i++)
		{
			FilterBankMel[i] = FilterBankMel[0] + (i * MelSpace);
			FilterBankHz[i] = GetFhzFromMel(FilterBankMel[i]);
		}

		//round to fft bins
		for (i = 0; i<(NumOfFilters + 2); i++)
			FilterBankFFTbins[i] = std::floor(((((double)FFT_FullSize + 1)*FilterBankHz[i]) / (double)SamplingRateHz));

		// fill the actual filter bank
		for (FilterIndex = 1; FilterIndex <=  NumOfFilters; FilterIndex++)//iterate over filters`
		{
			for (i = 0; i <FilterLength; i++)//iterate over filter points
			{
				if (i < FilterBankFFTbins[FilterIndex - 1])
				{
		 
					Output[MACRO_GetArrayIndex(FilterIndex - 1, i, FilterLength )] = 0;
				}
				else if ((FilterBankFFTbins[FilterIndex - 1] <= i) && (i <= FilterBankFFTbins[FilterIndex]))
				{
					Output[MACRO_GetArrayIndex(FilterIndex - 1, i, FilterLength)] = ((i - FilterBankFFTbins[FilterIndex - 1]) /\
						((FilterBankFFTbins[FilterIndex] - FilterBankFFTbins[FilterIndex - 1])));
				}
				else if ((FilterBankFFTbins[FilterIndex] <= i) && (i <= FilterBankFFTbins[FilterIndex + 1]))
				{
					Output[MACRO_GetArrayIndex(FilterIndex - 1, i, FilterLength)] = ((FilterBankFFTbins[FilterIndex + 1] - i) /\
						((FilterBankFFTbins[FilterIndex + 1] - FilterBankFFTbins[FilterIndex])));
				}
				else if (FilterBankFFTbins[FilterIndex + 1] < i)
				{
					Output[MACRO_GetArrayIndex(FilterIndex - 1, i, FilterLength)] = 0;
				}
			}//end of inner for loop
		}//end of outter for loop

		delete FilterBankMel;
		delete FilterBankHz;
		delete FilterBankFFTbins;

	}//end if

}// end func