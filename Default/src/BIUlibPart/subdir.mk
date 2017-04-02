################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/BIUlibPart/OptionParser.cpp 

OBJS += \
./src/BIUlibPart/OptionParser.o 

CPP_DEPS += \
./src/BIUlibPart/OptionParser.d 


# Each subdirectory must supply rules for building sources it contributes
src/BIUlibPart/%.o: ../src/BIUlibPart/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


