################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/post_proc.cu \

OBJS += \
./src/post_proc.o 

# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	nvcc -std=c++11 -Xcompiler -fopenmp -include ../../general_math.cu -include ../../format.cu -include ../../Space.cu -O3 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_50,code=compute_50 -ccbin g++ -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


