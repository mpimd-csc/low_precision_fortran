#include <stdio.h>
#include <stdint.h>

// Function to read the ID_AA64ISAR0_EL1 system register
static uint64_t read_id_aa64isar0_el1() {
    uint64_t value;
    asm volatile ("mrs %0, ID_AA64ISAR0_EL1" : "=r" (value));
    return value;
}

// Function to read the ID_AA64ISAR1_EL1 system register
static uint64_t read_id_aa64isar1_el1() {
    uint64_t value;
    asm volatile ("mrs %0, ID_AA64ISAR1_EL1" : "=r" (value));
    return value;
}

// Function to check for F16C flag (ARM equivalent is FP16)
static int cpu_has_fp16() {
    uint64_t id_aa64isar0 = read_id_aa64isar0_el1();
    // Bits 16-19 in ID_AA64ISAR0_EL1 indicate FP16 support
    uint64_t fp16_support = (id_aa64isar0 >> 16) & 0xF;
    return fp16_support != 0;
}
// Function to check for AVX512BF16 flag (ARM equivalent is BF16)
static int cpu_has_bf16() {
    uint64_t id_aa64isar1 = read_id_aa64isar1_el1();
    // Bits 4-7 in ID_AA64ISAR1_EL1 indicate BF16 support
    uint64_t bf16_support = (id_aa64isar1 >> 4) & 0xF;
    return bf16_support != 0;
}

int lpf_cpu_has_basic_fp16()
{
	return cpu_has_fp16();
}

int lpf_cpu_has_full_fp16() {
	return cpu_has_fp16();
}

int lpf_cpu_has_basic_bf16() {
	return cpu_has_bf16();
}

int lpf_cpu_has_full_bf16() {
	return cpu_has_bf16();
}


