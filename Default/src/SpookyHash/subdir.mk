################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/SpookyHash/SpookyV2.cpp 

OBJS += \
./src/SpookyHash/SpookyV2.o 

CPP_DEPS += \
./src/SpookyHash/SpookyV2.d 


# Each subdirectory must supply rules for building sources it contributes
src/SpookyHash/%.o: ../src/SpookyHash/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


