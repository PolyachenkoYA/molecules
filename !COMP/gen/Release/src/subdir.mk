################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/gen.cu \

OBJS += \
./src/gen.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-10.1/bin/nvcc -std=c++11 -Xcompiler -fopenmp -include ../../Space.cu -include ../../format.cu -include ../../general_math.cu -O3 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_50,code=compute_50 -ccbin g++ -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


