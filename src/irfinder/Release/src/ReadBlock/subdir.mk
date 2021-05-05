################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/ReadBlock/CoverageBlocks.cpp \
../src/ReadBlock/ReadBlockProcessor.cpp 

OBJS += \
./src/ReadBlock/CoverageBlocks.o \
./src/ReadBlock/ReadBlockProcessor.o 

CPP_DEPS += \
./src/ReadBlock/CoverageBlocks.d \
./src/ReadBlock/ReadBlockProcessor.d 


# Each subdirectory must supply rules for building sources it contributes
src/ReadBlock/%.o: ../src/ReadBlock/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


