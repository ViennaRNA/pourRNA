################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Flooder.cpp \
../src/GlobalParameter.cpp \
../src/MyState.cpp \
../src/ParKin_Explore.cpp \
../src/PriorityQueue.cpp \
../src/SC_DotPlot.cpp \
../src/SC_PartitionFunction.cpp \
../src/StateCollector.cpp \
../src/StatePairCollector.cpp \
../src/StructureUtils.cpp \
../src/TypeID.cpp \
../src/WalkGradientHashed.cpp \
../src/parKin.cpp 

OBJS += \
./src/Flooder.o \
./src/GlobalParameter.o \
./src/MyState.o \
./src/ParKin_Explore.o \
./src/PriorityQueue.o \
./src/SC_DotPlot.o \
./src/SC_PartitionFunction.o \
./src/StateCollector.o \
./src/StatePairCollector.o \
./src/StructureUtils.o \
./src/TypeID.o \
./src/WalkGradientHashed.o \
./src/parKin.o 

CPP_DEPS += \
./src/Flooder.d \
./src/GlobalParameter.d \
./src/MyState.d \
./src/ParKin_Explore.d \
./src/PriorityQueue.d \
./src/SC_DotPlot.d \
./src/SC_PartitionFunction.d \
./src/StateCollector.d \
./src/StatePairCollector.d \
./src/StructureUtils.d \
./src/TypeID.d \
./src/WalkGradientHashed.d \
./src/parKin.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


