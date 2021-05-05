################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Blocks/BAM2blocks.cpp \
../src/Blocks/CoverageBlock.cpp \
../src/Blocks/FragmentBlocks.cpp 

OBJS += \
./src/Blocks/BAM2blocks.o \
./src/Blocks/CoverageBlock.o \
./src/Blocks/FragmentBlocks.o 

CPP_DEPS += \
./src/Blocks/BAM2blocks.d \
./src/Blocks/CoverageBlock.d \
./src/Blocks/FragmentBlocks.d 


# Each subdirectory must supply rules for building sources it contributes
src/Blocks/%.o: ../src/Blocks/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


