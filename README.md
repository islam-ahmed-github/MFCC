# **MFCC**
Implementation of MFCC ( DSP algorithm: melody frequency cepstrum coefficients) features of a multi-dimentional data.\
it also support the calculation on DCT (discrete cosine transform), Acceleration, and rate of change.\
written in cpp with minimal modification required to use the methods in c.


# **DEPENDENCIES**
NONE.


# **INSTALLATION**
1- clone the repository.\
2- include the header in the project.



# **USAGE**
1- create an instance of the `MFCC` class.\
2- call the method `MFCC::GenerateFiltersBank` to generate the filter bank coefficients, to eleminate the most of multiplication operations in the runtime.\
3- call the function: `MFCC::GetMFCC` to get the mfcc features of the data. \
4- call the function: `MFCC::GetDCT` to get the discrete cosine transform of the data. \
5- call the function: `MFCC::GetDeltaDelta` to the get Acceleration (Delta), and the rate of change (Delta(Delta)).
