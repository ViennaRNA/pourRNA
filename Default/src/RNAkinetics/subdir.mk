################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/RNAkinetics/RNA_NeighMinFilter.cpp \
../src/RNAkinetics/RateMatrixUtil.cpp 

OBJS += \
./src/RNAkinetics/RNA_NeighMinFilter.o \
./src/RNAkinetics/RateMatrixUtil.o 

CPP_DEPS += \
./src/RNAkinetics/RNA_NeighMinFilter.d \
./src/RNAkinetics/RateMatrixUtil.d 


# Each subdirectory must supply rules for building sources it contributes
src/RNAkinetics/%.o: ../src/RNAkinetics/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


