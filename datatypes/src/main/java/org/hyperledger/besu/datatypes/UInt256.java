/*
 * Copyright ConsenSys AG.
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
 * an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations under the License.
 *
 * SPDX-License-Identifier: Apache-2.0
 */
package org.hyperledger.besu.datatypes;

public final class UInt256 {
  // 8 limbs, little-endian: x0 = least significant 32 bits
  private final int x0, x1, x2, x3, x4, x5, x6, x7;

  // ---- constructors ----

  public UInt256(final byte[] be) {
    if (be.length != 32) throw new IllegalArgumentException("Need 32 bytes");

    x7 = ((be[0] & 0xFF) << 24) | ((be[1] & 0xFF) << 16) | ((be[2] & 0xFF) << 8) | (be[3] & 0xFF);
    x6 = ((be[4] & 0xFF) << 24) | ((be[5] & 0xFF) << 16) | ((be[6] & 0xFF) << 8) | (be[7] & 0xFF);
    x5 = ((be[8] & 0xFF) << 24) | ((be[9] & 0xFF) << 16) | ((be[10] & 0xFF) << 8) | (be[11] & 0xFF);
    x4 =
        ((be[12] & 0xFF) << 24)
            | ((be[13] & 0xFF) << 16)
            | ((be[14] & 0xFF) << 8)
            | (be[15] & 0xFF);
    x3 =
        ((be[16] & 0xFF) << 24)
            | ((be[17] & 0xFF) << 16)
            | ((be[18] & 0xFF) << 8)
            | (be[19] & 0xFF);
    x2 =
        ((be[20] & 0xFF) << 24)
            | ((be[21] & 0xFF) << 16)
            | ((be[22] & 0xFF) << 8)
            | (be[23] & 0xFF);
    x1 =
        ((be[24] & 0xFF) << 24)
            | ((be[25] & 0xFF) << 16)
            | ((be[26] & 0xFF) << 8)
            | (be[27] & 0xFF);
    x0 =
        ((be[28] & 0xFF) << 24)
            | ((be[29] & 0xFF) << 16)
            | ((be[30] & 0xFF) << 8)
            | (be[31] & 0xFF);
  }

  public UInt256(
      final int x0,
      final int x1,
      final int x2,
      final int x3,
      final int x4,
      final int x5,
      final int x6,
      final int x7) {
    this.x0 = x0;
    this.x1 = x1;
    this.x2 = x2;
    this.x3 = x3;
    this.x4 = x4;
    this.x5 = x5;
    this.x6 = x6;
    this.x7 = x7;
  }

  public UInt256(final int x0) {
    this.x0 = x0;
    this.x1 = 0;
    this.x2 = 0;
    this.x3 = 0;
    this.x4 = 0;
    this.x5 = 0;
    this.x6 = 0;
    this.x7 = 0;
  }

  public UInt256(final long xl) {
    this.x0 = (int) xl;
    this.x1 = (int) (xl >>> 32);
    this.x2 = 0;
    this.x3 = 0;
    this.x4 = 0;
    this.x5 = 0;
    this.x6 = 0;
    this.x7 = 0;
  }

  public static UInt256 zero() {
    return new UInt256(0, 0, 0, 0, 0, 0, 0, 0);
  }

  // ---- conversion ----

  public int intValue() {
    return x0;
  }

  public long longValue() {
    return (x0 & 0xFFFFFFFFL) | ((x1 & 0xFFFFFFFFL) << 32);
  }

  public byte[] toBytesBE() {
    byte[] out = new byte[32];
    encodeInt(out, 0, x7);
    encodeInt(out, 4, x6);
    encodeInt(out, 8, x5);
    encodeInt(out, 12, x4);
    encodeInt(out, 16, x3);
    encodeInt(out, 20, x2);
    encodeInt(out, 24, x1);
    encodeInt(out, 28, x0);
    return out;
  }

  private static void encodeInt(final byte[] out, final int offset, final int v) {
    out[offset] = (byte) (v >>> 24);
    out[offset + 1] = (byte) (v >>> 16);
    out[offset + 2] = (byte) (v >>> 8);
    out[offset + 3] = (byte) v;
  }

  // ---- comparison ----

  private static int compare(final UInt256 a, final UInt256 b) {
    long[] aa = a.asLongs();
    long[] bb = b.asLongs();
    for (int i = 7; i >= 0; i--) {
      if (aa[i] < bb[i]) return -1;
      if (aa[i] > bb[i]) return 1;
    }
    return 0;
  }

  private long[] asLongs() {
    return new long[] {
      x0 & 0xFFFFFFFFL, x1 & 0xFFFFFFFFL,
      x2 & 0xFFFFFFFFL, x3 & 0xFFFFFFFFL,
      x4 & 0xFFFFFFFFL, x5 & 0xFFFFFFFFL,
      x6 & 0xFFFFFFFFL, x7 & 0xFFFFFFFFL
    };
  }

  private boolean isZero() {
    return (x0 | x1 | x2 | x3 | x4 | x5 | x6 | x7) == 0;
  }

  public int length() {
    if (x7 != 0) return 8;
    if (x6 != 0) return 7;
    if (x5 != 0) return 6;
    if (x4 != 0) return 5;
    if (x3 != 0) return 4;
    if (x2 != 0) return 3;
    if (x1 != 0) return 2;
    if (x0 != 0) return 1;
    return 0; // means isZero()
  }

  // ---- modulus ----
  // Simplified Knuth-style long division specialized for 256-bit ÷ 256-bit

  public UInt256 mod(final UInt256 b) {
    if (b.isZero()) throw new ArithmeticException("mod by zero");
    int cmp = compare(this, b);
    if (cmp < 0) return this;
    if (cmp == 0) return zero();

    // Estimate quotient from top limbs
    long aHi = x7 & 0xFFFFFFFFL;
    long aNext = x6 & 0xFFFFFFFFL;
    long bHi = b.x7 & 0xFFFFFFFFL;

    long dividend = (aHi << 32) | aNext;
    long qhat = dividend / bHi;
    if (qhat > 0xFFFFFFFFL) qhat = 0xFFFFFFFFL;

    // Copy remainder
    long r0 = x0 & 0xFFFFFFFFL, r1 = x1 & 0xFFFFFFFFL, r2 = x2 & 0xFFFFFFFFL, r3 = x3 & 0xFFFFFFFFL;
    long r4 = x4 & 0xFFFFFFFFL, r5 = x5 & 0xFFFFFFFFL, r6 = x6 & 0xFFFFFFFFL, r7 = x7 & 0xFFFFFFFFL;

    // Multiply-subtract remainder -= qhat * divisor
    long carry = 0, borrow = 0;

    long prod = (b.x0 & 0xFFFFFFFFL) * qhat + carry;
    carry = prod >>> 32;
    long diff = r0 - (prod & 0xFFFFFFFFL) - borrow;
    r0 = diff & 0xFFFFFFFFL;
    borrow = (diff >> 63) & 1;

    prod = (b.x1 & 0xFFFFFFFFL) * qhat + carry;
    carry = prod >>> 32;
    diff = r1 - (prod & 0xFFFFFFFFL) - borrow;
    r1 = diff & 0xFFFFFFFFL;
    borrow = (diff >> 63) & 1;

    prod = (b.x2 & 0xFFFFFFFFL) * qhat + carry;
    carry = prod >>> 32;
    diff = r2 - (prod & 0xFFFFFFFFL) - borrow;
    r2 = diff & 0xFFFFFFFFL;
    borrow = (diff >> 63) & 1;

    prod = (b.x3 & 0xFFFFFFFFL) * qhat + carry;
    carry = prod >>> 32;
    diff = r3 - (prod & 0xFFFFFFFFL) - borrow;
    r3 = diff & 0xFFFFFFFFL;
    borrow = (diff >> 63) & 1;

    prod = (b.x4 & 0xFFFFFFFFL) * qhat + carry;
    carry = prod >>> 32;
    diff = r4 - (prod & 0xFFFFFFFFL) - borrow;
    r4 = diff & 0xFFFFFFFFL;
    borrow = (diff >> 63) & 1;

    prod = (b.x5 & 0xFFFFFFFFL) * qhat + carry;
    carry = prod >>> 32;
    diff = r5 - (prod & 0xFFFFFFFFL) - borrow;
    r5 = diff & 0xFFFFFFFFL;
    borrow = (diff >> 63) & 1;

    prod = (b.x6 & 0xFFFFFFFFL) * qhat + carry;
    carry = prod >>> 32;
    diff = r6 - (prod & 0xFFFFFFFFL) - borrow;
    r6 = diff & 0xFFFFFFFFL;
    borrow = (diff >> 63) & 1;

    prod = (b.x7 & 0xFFFFFFFFL) * qhat + carry;
    // carry = prod >>> 32;
    diff = r7 - (prod & 0xFFFFFFFFL) - borrow;
    r7 = diff & 0xFFFFFFFFL;
    // borrow = (diff >> 63) & 1;

    UInt256 rem =
        new UInt256((int) r0, (int) r1, (int) r2, (int) r3, (int) r4, (int) r5, (int) r6, (int) r7);

    // Ensure remainder < divisor
    while (compare(rem, b) >= 0) {
      rem = rem.sub(b);
    }
    return rem;
  }

  // ---- subtraction (this - b), assuming this >= b ----

  private UInt256 sub(final UInt256 b) {
    long r0 = (x0 & 0xFFFFFFFFL) - (b.x0 & 0xFFFFFFFFL);
    long borrow = (r0 >> 63) & 1;
    r0 &= 0xFFFFFFFFL;

    long r1 = (x1 & 0xFFFFFFFFL) - (b.x1 & 0xFFFFFFFFL) - borrow;
    borrow = (r1 >> 63) & 1;
    r1 &= 0xFFFFFFFFL;

    long r2 = (x2 & 0xFFFFFFFFL) - (b.x2 & 0xFFFFFFFFL) - borrow;
    borrow = (r2 >> 63) & 1;
    r2 &= 0xFFFFFFFFL;

    long r3 = (x3 & 0xFFFFFFFFL) - (b.x3 & 0xFFFFFFFFL) - borrow;
    borrow = (r3 >> 63) & 1;
    r3 &= 0xFFFFFFFFL;

    long r4 = (x4 & 0xFFFFFFFFL) - (b.x4 & 0xFFFFFFFFL) - borrow;
    borrow = (r4 >> 63) & 1;
    r4 &= 0xFFFFFFFFL;

    long r5 = (x5 & 0xFFFFFFFFL) - (b.x5 & 0xFFFFFFFFL) - borrow;
    borrow = (r5 >> 63) & 1;
    r5 &= 0xFFFFFFFFL;

    long r6 = (x6 & 0xFFFFFFFFL) - (b.x6 & 0xFFFFFFFFL) - borrow;
    borrow = (r6 >> 63) & 1;
    r6 &= 0xFFFFFFFFL;

    long r7 = (x7 & 0xFFFFFFFFL) - (b.x7 & 0xFFFFFFFFL) - borrow;
    r7 &= 0xFFFFFFFFL;

    return new UInt256(
        (int) r0, (int) r1, (int) r2, (int) r3, (int) r4, (int) r5, (int) r6, (int) r7);
  }

  // ---- debugging ----

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder("0x");
    for (byte b : toBytesBE()) {
      sb.append(String.format("%02x", b));
    }
    return sb.toString();
  }
}
