

#pragma once

 

#include<stdint.h>
#include<cmath>
#include<windows.h>
#include<limits.h>
#include<float.h>


 
#define MACRO_PI (3.14159265358979323846)
#define MFCCMACRO_FilterWidth  (((MFCCMACRO_FFT_FullSize)/2)+1)
#define MACRO_GetArrayIndex(RowIndex,ColIndex,RowLength)  (((RowIndex)*RowLength)+(ColIndex))

class MFCC
{

public:

	MFCC();

	virtual ~MFCC();

	void MFCC::GenerateFiltersBank(double *Output, UINT32 NumOfFilters, UINT32 SamplingRateHz, UINT32 UpperFreqHz, UINT32 LowerFreqHz, UINT32 FFT_FullSize);
	void MFCC::GetMFCC(const double *InputFrame, const uint32_t InputFrameLength, const double *FiltersBank, UINT32 NumOfFilters, double *OutputFrame);
	void MFCC::GetDCT(const double *InputFrame, uint32_t InputFrameLength, double *OutputFrame, uint32_t NumOfCoeff);
	void MFCC::GetDeltaDelta(const double *InputMfcc, uint32_t NumOfVectors, uint32_t VectorLength, double *OutputMfccDD);
	 
	void MFCC::GetMFCCDCT(const double *InputFrame, const uint32_t InputFrameLength, const double *FiltersBank, UINT32 NumOfFilters, UINT32 NumOfCoefficients, double *OutputFrame);


protected:



private:


};

 

