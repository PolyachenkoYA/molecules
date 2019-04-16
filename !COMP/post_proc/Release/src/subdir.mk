################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/post_proc.cu 

OBJS += \
./src/post_proc.o 

CU_DEPS += \
./src/post_proc.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-9.2/bin/nvcc -include ../../Space.cu -include ../../general_math.cu -include ../../format.cu -O3 -Xcompiler -fopenmp -std=c++11 -gencode arch=compute_50,code=sm_50  -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-9.2/bin/nvcc -include ../../Space.cu -include ../../general_math.cu -include ../../format.cu -O3 -Xcompiler -fopenmp -std=c++11 --compile --relocatable-device-code=false -gencode arch=compute_50,code=compute_50 -gencode arch=compute_50,code=sm_50  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


