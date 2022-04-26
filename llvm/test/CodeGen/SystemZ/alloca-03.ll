; NOTE: Assertions have been autogenerated by utils/update_llc_test_checks.py
; RUN: llc < %s -mtriple=s390x-linux-gnu | FileCheck %s

; Allocate 8 bytes, no need to align stack.
define void @f0() {
; CHECK-LABEL: f0:
; CHECK:       # %bb.0:
; CHECK-NEXT:    aghi %r15, -168
; CHECK-NEXT:    .cfi_def_cfa_offset 328
; CHECK-NEXT:    mvghi 160(%r15), 10
; CHECK-NEXT:    aghi %r15, 168
; CHECK-NEXT:    br %r14
  %x = alloca i64
  store volatile i64 10, i64* %x
  ret void
}

; Allocate %len * 8, no need to align stack.
define void @f1(i64 %len) {
; CHECK-LABEL: f1:
; CHECK:       # %bb.0:
; CHECK-NEXT:    stmg %r11, %r15, 88(%r15)
; CHECK-NEXT:    .cfi_offset %r11, -72
; CHECK-NEXT:    .cfi_offset %r15, -40
; CHECK-NEXT:    aghi %r15, -160
; CHECK-NEXT:    .cfi_def_cfa_offset 320
; CHECK-NEXT:    lgr %r11, %r15
; CHECK-NEXT:    .cfi_def_cfa_register %r11
; CHECK-NEXT:    lgr %r1, %r15
; CHECK-NEXT:    sllg %r0, %r2, 3
; CHECK-NEXT:    sgr %r1, %r0
; CHECK-NEXT:    la %r2, 160(%r1)
; CHECK-NEXT:    lgr %r15, %r1
; CHECK-NEXT:    mvghi 0(%r2), 10
; CHECK-NEXT:    lmg %r11, %r15, 248(%r11)
; CHECK-NEXT:    br %r14
  %x = alloca i64, i64 %len
  store volatile i64 10, i64* %x
  ret void
}

; Static alloca, align 128.
define void @f2() {
; CHECK-LABEL: f2:
; CHECK:       # %bb.0:
; CHECK-NEXT:    stmg %r11, %r15, 88(%r15)
; CHECK-NEXT:    .cfi_offset %r11, -72
; CHECK-NEXT:    .cfi_offset %r15, -40
; CHECK-NEXT:    aghi %r15, -160
; CHECK-NEXT:    .cfi_def_cfa_offset 320
; CHECK-NEXT:    lgr %r11, %r15
; CHECK-NEXT:    .cfi_def_cfa_register %r11
; CHECK-NEXT:    lgr %r1, %r15
; CHECK-NEXT:    aghi %r1, -128
; CHECK-NEXT:    la %r2, 280(%r1)
; CHECK-NEXT:    nill %r2, 65408
; CHECK-NEXT:    lgr %r15, %r1
; CHECK-NEXT:    mvghi 0(%r2), 10
; CHECK-NEXT:    lmg %r11, %r15, 248(%r11)
; CHECK-NEXT:    br %r14
  %x = alloca i64, i64 1, align 128
  store volatile i64 10, i64* %x, align 128
  ret void
}

; Dynamic alloca, align 128.
define void @f3(i64 %len) {
; CHECK-LABEL: f3:
; CHECK:       # %bb.0:
; CHECK-NEXT:    stmg %r11, %r15, 88(%r15)
; CHECK-NEXT:    .cfi_offset %r11, -72
; CHECK-NEXT:    .cfi_offset %r15, -40
; CHECK-NEXT:    aghi %r15, -160
; CHECK-NEXT:    .cfi_def_cfa_offset 320
; CHECK-NEXT:    lgr %r11, %r15
; CHECK-NEXT:    .cfi_def_cfa_register %r11
; CHECK-NEXT:    lgr %r1, %r15
; CHECK-NEXT:    sllg %r0, %r2, 3
; CHECK-NEXT:    sgr %r1, %r0
; CHECK-NEXT:    lay %r15, -120(%r1)
; CHECK-NEXT:    la %r1, 160(%r1)
; CHECK-NEXT:    nill %r1, 65408
; CHECK-NEXT:    mvghi 0(%r1), 10
; CHECK-NEXT:    lmg %r11, %r15, 248(%r11)
; CHECK-NEXT:    br %r14
  %x = alloca i64, i64 %len, align 128
  store volatile i64 10, i64* %x, align 128
  ret void
}

; Static alloca w/out alignment - part of frame.
define void @f4() {
; CHECK-LABEL: f4:
; CHECK:       # %bb.0:
; CHECK-NEXT:    aghi %r15, -168
; CHECK-NEXT:    .cfi_def_cfa_offset 328
; CHECK-NEXT:    mvhi 164(%r15), 10
; CHECK-NEXT:    aghi %r15, 168
; CHECK-NEXT:    br %r14
  %x = alloca i32
  store volatile i32 10, i32* %x
  ret void
}

; Static alloca of one i32, aligned by 128.
define void @f5() {
; CHECK-LABEL: f5:
; CHECK:       # %bb.0:
; CHECK-NEXT:    stmg %r11, %r15, 88(%r15)
; CHECK-NEXT:    .cfi_offset %r11, -72
; CHECK-NEXT:    .cfi_offset %r15, -40
; CHECK-NEXT:    aghi %r15, -160
; CHECK-NEXT:    .cfi_def_cfa_offset 320
; CHECK-NEXT:    lgr %r11, %r15
; CHECK-NEXT:    .cfi_def_cfa_register %r11
; CHECK-NEXT:    lgr %r1, %r15
; CHECK-NEXT:    aghi %r1, -128
; CHECK-NEXT:    la %r2, 280(%r1)
; CHECK-NEXT:    nill %r2, 65408
; CHECK-NEXT:    lgr %r15, %r1
; CHECK-NEXT:    mvhi 0(%r2), 10
; CHECK-NEXT:    lmg %r11, %r15, 248(%r11)
; CHECK-NEXT:    br %r14
  %x = alloca i32, i64 1, align 128
  store volatile i32 10, i32* %x
  ret void
}
