################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/gen.cu 

OBJS += \
./src/gen.o 

CU_DEPS += \
./src/gen.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-9.2/bin/nvcc -include ../../general_math.cu -include ../../format.cu -include ../../Space.cu -G -g -O0 -Xcompiler -fopenmp -std=c++11   -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-9.2/bin/nvcc -include ../../general_math.cu -include ../../format.cu -include ../../Space.cu -G -g -O0 -Xcompiler -fopenmp -std=c++11 --compile --relocatable-device-code=false  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

