################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -f

# All of the sources participating in the build are defined here
-include sources.mk
-include src/subdir.mk
-include subdir.mk
-include objects.mk

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: post_proc

# Tool invocations
post_proc: $(OBJS) $(USER_OBJS) $(DEPS)
	@echo 'Building target: $@'
	@echo 'Invoking: NVCC linker'
	nvcc --cudart=static -ccbin g++ -gencode arch=compute_50,code=sm_50 -gencode arch=compute_50,code=compute_50 -o "post_proc" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(EXECUTABLES)$(OBJS) post_proc
	-@echo ' '

.PHONY: all clean dependents

-include ../makefile.targets
