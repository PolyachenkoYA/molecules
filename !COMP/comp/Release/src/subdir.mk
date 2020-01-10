################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/comp.cu 

OBJS += \
./src/comp.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-10.1/bin/nvcc -std=c++11 -Xcompiler -fopenmp -include ../../general_math.cu -include ../../format.cu -include ../../Space.cu -O3 --use_fast_math -gencode arch=compute_50,code=sm_50 -gencode arch=compute_50,code=compute_50 -ccbin g++ -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


