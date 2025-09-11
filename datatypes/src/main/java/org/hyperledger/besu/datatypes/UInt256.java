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

import java.util.Arrays;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.LongVector;
import jdk.incubator.vector.VectorSpecies;

/**
 * 256-bits wide unsigned integer class.
 *
 * <p>This class is an optimised version of BigInteger for fixed width 8byte integers
 */
public final class UInt256 {
  // Wraps a view over limbs array from 0..length.
  // Limbs are little-endian and can perhaps have more elements than length, which are ignored.
  private final int[] limbs;
  private final int length;

  // Constants
  private static final long BASE = 1L << 32;
  private static final long MASK32 = 0xFFFF_FFFFL;

  private static final VectorSpecies<Long> LSPEC = LongVector.SPECIES_PREFERRED;
  private static final VectorSpecies<Integer> ISPEC = IntVector.SPECIES_PREFERRED;

  // reciprocal32 lookup table, fits in 15bits, but last bit is always set.
  private static final short[] REC_LUT = new short[] {
    32737, 32673, 32609, 32546, 32483, 32420, 32357, 32295, 32233, 32171, 32109, 32048, 31987, 31926, 31865, 31805,
    31744, 31684, 31625, 31565, 31506, 31447, 31388, 31329, 31271, 31212, 31154, 31097, 31039, 30982, 30924, 30868,
    30811, 30754, 30698, 30642, 30586, 30530, 30475, 30419, 30364, 30309, 30255, 30200, 30146, 30092, 30038, 29984,
    29930, 29877, 29824, 29771, 29718, 29666, 29613, 29561, 29509, 29457, 29405, 29354, 29303, 29251, 29200, 29150,
    29099, 29049, 28998, 28948, 28898, 28849, 28799, 28750, 28700, 28651, 28602, 28554, 28505, 28457, 28409, 28360,
    28313, 28265, 28217, 28170, 28123, 28075, 28029, 27982, 27935, 27889, 27842, 27796, 27750, 27704, 27658, 27613,
    27568, 27522, 27477, 27432, 27387, 27343, 27298, 27254, 27209, 27165, 27121, 27078, 27034, 26990, 26947, 26904,
    26861, 26818, 26775, 26732, 26690, 26647, 26605, 26563, 26521, 26479, 26437, 26395, 26354, 26312, 26271, 26230,
    26189, 26148, 26108, 26067, 26026, 25986, 25946, 25906, 25866, 25826, 25786, 25747, 25707, 25668, 25628, 25589,
    25550, 25511, 25473, 25434, 25395, 25357, 25319, 25281, 25242, 25205, 25167, 25129, 25091, 25054, 25016, 24979,
    24942, 24905, 24868, 24831, 24794, 24758, 24721, 24685, 24649, 24612, 24576, 24540, 24504, 24469, 24433, 24397,
    24362, 24327, 24291, 24256, 24221, 24186, 24151, 24117, 24082, 24047, 24013, 23979, 23944, 23910, 23876, 23842,
    23808, 23774, 23741, 23707, 23674, 23640, 23607, 23574, 23541, 23508, 23475, 23442, 23409, 23377, 23344, 23312,
    23279, 23247, 23215, 23183, 23151, 23119, 23087, 23055, 23023, 22992, 22960, 22929, 22898, 22866, 22835, 22804,
    22773, 22742, 22711, 22681, 22650, 22619, 22589, 22559, 22528, 22498, 22468, 22438, 22408, 22378, 22348, 22318,
    22289, 22259, 22229, 22200, 22171, 22141, 22112, 22083, 22054, 22025, 21996, 21967, 21938, 21910, 21881, 21853,
    21824, 21796, 21767, 21739, 21711, 21683, 21655, 21627, 21599, 21571, 21544, 21516, 21488, 21461, 21433, 21406,
    21379, 21352, 21324, 21297, 21270, 21243, 21216, 21190, 21163, 21136, 21110, 21083, 21056, 21030, 21004, 20977,
    20951, 20925, 20899, 20873, 20847, 20821, 20795, 20769, 20744, 20718, 20693, 20667, 20642, 20616, 20591, 20566,
    20540, 20515, 20490, 20465, 20440, 20415, 20390, 20366, 20341, 20316, 20292, 20267, 20243, 20218, 20194, 20170,
    20145, 20121, 20097, 20073, 20049, 20025, 20001, 19977, 19953, 19930, 19906, 19882, 19859, 19835, 19812, 19789,
    19765, 19742, 19719, 19696, 19672, 19649, 19626, 19603, 19581, 19558, 19535, 19512, 19489, 19467, 19444, 19422,
    19399, 19377, 19354, 19332, 19310, 19288, 19265, 19243, 19221, 19199, 19177, 19155, 19133, 19112, 19090, 19068,
    19046, 19025, 19003, 18982, 18960, 18939, 18917, 18896, 18875, 18854, 18832, 18811, 18790, 18769, 18748, 18727,
    18706, 18686, 18665, 18644, 18623, 18603, 18582, 18561, 18541, 18520, 18500, 18479, 18459, 18439, 18419, 18398,
    18378, 18358, 18338, 18318, 18298, 18278, 18258, 18238, 18218, 18199, 18179, 18159, 18139, 18120, 18100, 18081,
    18061, 18042, 18022, 18003, 17984, 17964, 17945, 17926, 17907, 17888, 17869, 17850, 17831, 17812, 17793, 17774,
    17755, 17736, 17718, 17699, 17680, 17662, 17643, 17624, 17606, 17587, 17569, 17551, 17532, 17514, 17496, 17477,
    17459, 17441, 17423, 17405, 17387, 17369, 17351, 17333, 17315, 17297, 17279, 17261, 17244, 17226, 17208, 17191,
    17173, 17155, 17138, 17120, 17103, 17085, 17068, 17051, 17033, 17016, 16999, 16982, 16964, 16947, 16930, 16913,
    16896, 16879, 16862, 16845, 16828, 16811, 16794, 16778, 16761, 16744, 16727, 16711, 16694, 16677, 16661, 16644,
    16628, 16611, 16595, 16578, 16562, 16546, 16529, 16513, 16497, 16481, 16464, 16448, 16432, 16416, 16400, 16384 };

  // --- Getters for testing ---
  int length() {
    return length;
  }

  int[] limbs() {
    return limbs;
  }

  // --- Preallocating small integers 0..nSmallInts ---
  private static final int nSmallInts = 256;
  private static final UInt256[] smallInts = new UInt256[nSmallInts];

  static {
    smallInts[0] = new UInt256(new int[] {});
    for (int i = 1; i < nSmallInts; i++) {
      smallInts[i] = new UInt256(new int[] {i});
    }
  }

  /** The constant 0. */
  public static final UInt256 ZERO = smallInts[0];

  /** The constant 1. */
  public static final UInt256 ONE = smallInts[1];

  /** The constant 2. */
  public static final UInt256 TWO = smallInts[2];

  /** The constant 10. */
  public static final UInt256 TEN = smallInts[10];

  /** The constant 16. */
  public static final UInt256 SIXTEEN = smallInts[16];

  // --- Constructors ---

  UInt256(final int[] limbs, final int length) {
    this.limbs = limbs;
    this.length = length;
    // Unchecked length: assumes length is properly set.
  }

  UInt256(final int[] limbs) {
    int i = Math.min(limbs.length - 1, 7);
    while ((i >= 0) && (limbs[i] == 0)) i--;
    this.limbs = limbs;
    this.length = i + 1;
  }

  /**
   * Instantiates a new UInt256 from byte[].
   *
   * @param bytes raw bytes in BigEndian order.
   * @return Big-endian UInt256 represented by the bytes.
   */
  public static UInt256 fromBytesBE(final byte[] bytes) {
    int offset = 0;
    while ((offset < bytes.length) && (bytes[offset] == 0x00)) ++offset;
    int nBytes = bytes.length - offset;
    if (nBytes == 0) return ZERO;
    int len = (nBytes + 3) / 4;
    int[] limbs = new int[len];

    int i;
    int base;
    // Up to most significant limb take 4 bytes.
    for (i = 0, base = bytes.length - 4; i < len - 1; ++i, base = base - 4) {
      limbs[i] =
          (bytes[base] << 24)
              | ((bytes[base + 1] & 0xFF) << 16)
              | ((bytes[base + 2] & 0xFF) << 8)
              | ((bytes[base + 3] & 0xFF));
    }
    // Last effective limb
    limbs[i] =
        switch (nBytes - i * 4) {
          case 1 -> ((bytes[offset] & 0xFF));
          case 2 -> (((bytes[offset] & 0xFF) << 8) | (bytes[offset + 1] & 0xFF));
          case 3 ->
              (((bytes[offset] & 0xFF) << 16)
                  | ((bytes[offset + 1] & 0xFF) << 8)
                  | (bytes[offset + 2] & 0xFF));
          case 4 ->
              ((bytes[offset] << 24)
                  | ((bytes[offset + 1] & 0xFF) << 16)
                  | ((bytes[offset + 2] & 0xFF) << 8)
                  | (bytes[offset + 3] & 0xFF));
          default -> throw new IllegalStateException("Unexpected value");
        };
    return new UInt256(limbs, len);
  }

  /**
   * Instantiates a new UInt256 from an int.
   *
   * @param value int value to convert to UInt256.
   * @return The UInt256 equivalent of value.
   */
  public static UInt256 fromInt(final int value) {
    if (0 <= value && value < nSmallInts) return smallInts[value];
    return new UInt256(new int[] {value}, 1);
  }

  /**
   * Instantiates a new UInt256 from a long.
   *
   * @param value long value to convert to UInt256.
   * @return The UInt256 equivalent of value.
   */
  public static UInt256 fromLong(final long value) {
    if (0 <= value && value < nSmallInts) return smallInts[(int) value];
    return new UInt256(new int[] {(int) value, (int) (value >>> 32)});
  }

  // ---- conversion ----

  /**
   * Convert to int.
   *
   * @return Value truncated to an int, possibly lossy.
   */
  public int intValue() {
    return (limbs.length == 0 ? 0 : limbs[0]);
  }

  /**
   * Convert to long.
   *
   * @return Value truncated to a long, possibly lossy.
   */
  public long longValue() {
    switch (length) {
      case 0 -> {
        return 0L;
      }
      case 1 -> {
        return (limbs[0] & 0xFFFFFFFFL);
      }
      default -> {
        return (limbs[0] & 0xFFFFFFFFL) | ((limbs[1] & 0xFFFFFFFFL) << 32);
      }
    }
  }

  /**
   * Convert to BigEndian byte array.
   *
   * @return Big-endian ordered bytes for this UInt256 value.
   */
  public byte[] toBytesBE() {
    byte[] out = new byte[length * 4];
    for (int i = 0; i < length; ++i) {
      int offset = i * 4;
      int v = limbs[length - i - 1];
      out[offset] = (byte) (v >>> 24);
      out[offset + 1] = (byte) (v >>> 16);
      out[offset + 2] = (byte) (v >>> 8);
      out[offset + 3] = (byte) v;
    }
    return out;
  }

  // ---- comparison ----

  /**
   * Is the value 0 ?
   *
   * @return true if this UInt256 value is 0.
   */
  public boolean isZero() {
    if (length == 0) return true;
    for (int i = 0; i < length; i++) {
      if (limbs[i] != 0) return false;
    }
    return true;
  }

  /**
   * Compares two UInt256.
   *
   * @param a left UInt256
   * @param b right UInt256
   * @return 0 if a == b, negative if a &lt; b and positive if a &gt; b.
   */
  public static int compare(final UInt256 a, final UInt256 b) {
    int comp = Integer.compare(a.length, b.length);
    if (comp != 0) return comp;
    for (int i = a.length - 1; i >= 0; i--) {
      comp = Integer.compareUnsigned(a.limbs[i], b.limbs[i]);
      if (comp != 0) return comp;
    }
    return 0;
  }

  private int[] shiftLeftWithCarry(final int shift) {
    int limbShift = shift / 32;
    int bitShift = shift % 32;
    if (bitShift == 0) {
      return Arrays.copyOfRange(limbs, limbShift, length + limbShift);
    }

    int[] res = new int[length + limbShift + 1];
    int j = limbShift;
    int carry = 0;
    for (int i = 0; (i < length) && (j < 8); ++i, ++j) {
      res[j] = (limbs[i] << bitShift) | carry;
      carry = limbs[i] >>> (32 - bitShift);
    }
    res[j] = carry; // last carry
    return res;
  }

  /**
   * Shifts value to the left.
   *
   * @param shift number of bits to shift. If negative, shift right instead.
   * @return Shifted UInt256 value.
   */
  public UInt256 shiftLeft(final int shift) {
    if (shift < 0) return shiftRight(-shift);
    if (shift == 0 || isZero()) return this;
    if (shift >= length * 32) return ZERO;
    return new UInt256(shiftLeftWithCarry(shift));
  }

  /**
   * Shifts value to the right.
   *
   * @param shift number of bits to shift. If negative, shift left instead.
   * @return Shifted UInt256 value.
   */
  public UInt256 shiftRight(final int shift) {
    if (shift < 0) return shiftLeft(-shift);
    if (shift == 0 || isZero()) return this;
    if (shift >= length * 32) return ZERO;

    int limbShift = shift / 32;
    int bitShift = shift % 32;
    if (bitShift == 0) {
      return new UInt256(Arrays.copyOfRange(limbs, limbShift, length), length - limbShift);
    }

    int[] res = new int[this.limbs.length - limbShift];
    int j = this.limbs.length - 1 - limbShift; // res index
    int carry = 0;
    for (int i = this.limbs.length - 1; j >= 0; i--, j--) {
      int r = (limbs[i] >>> bitShift) | carry;
      res[j] = r;
      carry = limbs[i] << (32 - bitShift);
    }
    return new UInt256(res);
  }

  /*
  private int int_reciprocal(final int d) {
    // Assert that d < 0

    int d0 = d & 1;
    int d10 = d >>> 22;  // 10 msb
    int d21 = 1 + (d >> 11);
    int d31 = (d + 1) >> 1;  // Should be < MAX ?
    int v0 = REC_LUT[d10 - 512];
    
    long v0sq = (long) v0 * v0;
    int v0sqd21 = (int) ((v0sq * d21) >>> 32);
    int v1 = (v0 << 4) - v0sqd21 - 1;


    // int v2 = (v1 << 13) + ((v1 * ((1 << 60) - v1 * d40)) >> 47);
    // let e = ((v2 >> 1) & (ZERO - d0)) - v2 * d63;
    // let v3 = (mul_hi(v2, e) >> 1) + (v2 << 31);
    // let v4 = v3 - muladd_hi(v3, d, d) - d;

    // v4.0
  }
  */

  /**
   * Reduce modulo divisor.
   *
   * @param divisor The modulus of the reduction
   * @return The remainder modulo {@code divisor}.
   */
  public UInt256 modClassic(final UInt256 divisor) {
    if (divisor.isZero()) return ZERO;
    int cmp = compare(this, divisor);
    if (cmp < 0) return this;
    if (cmp == 0) return ZERO;

    int n = divisor.length;
    if (n == 1) {
      // divNby1Wide
      long d = divisor.limbs[0] & 0xFFFFFFFFL;
      long rem = 0;
      // Process from most significant limb downwards
      for (int i = length - 1; i >= 0; i--) {
        long cur = (rem << 32) | (limbs[i] & 0xFFFFFFFFL);
        rem = Long.remainderUnsigned(cur, d);
      }
      return new UInt256(new int[] {(int) rem});
    }

    // Case 64 bits ?

    // --- Knuth Division ---

    // Normalize
    // Makes 2 copies of limbs: optim ?
    int shift = Integer.numberOfLeadingZeros(divisor.limbs[n - 1]);
    int[] vLimbs = divisor.shiftLeftWithCarry(shift);
    int[] uLimbs = this.shiftLeftWithCarry(shift);

    // Main division loop
    int m = uLimbs.length - n - 1;
    long[] vLimbsWide = new long[vLimbs.length];
    for (int i = 0; i < vLimbs.length; i++) {
      vLimbsWide[i] = vLimbs[i] & 0xFFFFFFFFL;
    }

    long vn1 = vLimbsWide[n - 1];
    long vn2 = (n > 1) ? vLimbsWide[n - 2] : 0L;
    for (int j = m; j >= 0; j--) {
      // Estimate quotient
      long qhat = trialQuotient(uLimbs, vn1, vn2, j, n);
      int borrow = mulSub(uLimbs, j, vLimbsWide, n, qhat);
      if (borrow != 0) {
        addBack(uLimbs, j, vLimbsWide, n, 1);
        qhat -= 1;
      }
    }

    // Unnormalize remainder
    return (new UInt256(uLimbs, n)).shiftRight(shift);
  }

  /**
   * Reduce modulo divisor.
   *
   * @param divisor The modulus of the reduction
   * @return The remainder modulo {@code divisor}.
   */
  public UInt256 mod(final UInt256 divisor) {
    if (divisor.isZero()) return ZERO;
    int cmp = compare(this, divisor);
    if (cmp < 0) return this;
    if (cmp == 0) return ZERO;

    int n = divisor.length;
    if (n == 1) {
      // divNby1Wide
      long d = divisor.limbs[0] & 0xFFFFFFFFL;
      long rem = 0;
      // Process from most significant limb downwards
      for (int i = length - 1; i >= 0; i--) {
        long cur = (rem << 32) | (limbs[i] & 0xFFFFFFFFL);
        rem = Long.remainderUnsigned(cur, d);
      }
      return new UInt256(new int[] {(int) rem});
    }

    // Case 64 bits ?

    // --- Knuth Division ---

    // Normalize
    // Makes 2 copies of limbs: optim ?
    int shift = Integer.numberOfLeadingZeros(divisor.limbs[n - 1]);
    int[] vLimbs = divisor.shiftLeftWithCarry(shift);
    int[] uLimbs = this.shiftLeftWithCarry(shift);

    // Main division loop
    int m = uLimbs.length - n - 1;
    long[] vLimbsWide = new long[vLimbs.length];
    for (int i = 0; i < vLimbs.length; i++) {
      vLimbsWide[i] = vLimbs[i] & 0xFFFFFFFFL;
    }

    long vn1 = vLimbsWide[n - 1];
    long vn2 = (n > 1) ? vLimbsWide[n - 2] : 0L;
    for (int j = m; j >= 0; j--) {
      // Estimate quotient
      long qhat = trialQuotient(uLimbs, vn1, vn2, j, n);
      int borrow = mulSubSIMD(uLimbs, j, vLimbsWide, n, qhat);
      if (borrow != 0) {
        addBackSIMD(uLimbs, j, vLimbs, n, 1);
        qhat -= 1;
      }
    }

    // Unnormalize remainder
    return (new UInt256(uLimbs, n)).shiftRight(shift);
  }

  private static long trialQuotient(final int[] U, final long vn1, final long vn2, final int j, final int n) {
    long ujn = U[j + n] & MASK32;
    long ujn1 = U[j + n - 1] & MASK32;
    long ujn2 = U[j + n - 2] & MASK32;

    // div2by1 using widening: (U[j+n] * B + U[j+n-1]) / V[n-1]
    long numerator = (ujn << 32) | ujn1;
    long qhat = numerator / vn1;
    long rhat = numerator % vn1;

    if (qhat >= BASE) {
      qhat = BASE - 1;
      rhat = numerator - qhat * vn1; // recompute remainder
    }

    // Correction: check if qhat*vn2 > rhat*b + U[j+n-2]
    while (qhat * vn2 > (rhat << 32) + ujn2) {
      qhat -= 1;
      rhat += vn1;
      if (rhat >= BASE) break; // can't loop forever
    }

    return qhat;
  }

  private static int mulSub(final int[] U, final int j, final long[] VL, final int n, final long qhat) {
    long borrow = 0;
    for (int i = 0; i < n; i++) {
      long prod = VL[i] * qhat;
      long sub = (U[i + j] & 0xFFFFFFFFL) - (prod & 0xFFFFFFFFL) - borrow;
      U[i + j] = (int) sub;
      borrow = (prod >>> 32) - (sub >> 63);
    }
    long sub = (U[j + n] & 0xFFFFFFFFL) - borrow;
    U[j + n] = (int) sub;
    return (sub < 0) ? 1 : 0;
  }

  private static int mulSubSIMD(final int[] U, final int j, final long[] VL, final int n, final long qhat) {
    final int lanes = LSPEC.length();

    final long[] vBuf = new long[lanes];
    final long[] uBuf = new long[lanes];
    final long[] prodBuf = new long[lanes];

    long carry = 0;
    int i = 0;

    for (; i + lanes <= n; i += lanes) {
      for (int k = 0; k < lanes; k++) {
        vBuf[k] = VL[i + k];
        uBuf[k] = U[j + i + k] & MASK32;
      }
      // Use offset for vVec
      LongVector vVec = LongVector.fromArray(LSPEC, vBuf, 0);
      LongVector pVec = vVec.mul(qhat);
      pVec.intoArray(prodBuf, 0);

      for (int k = 0; k < lanes; k++) {
        long prod = prodBuf[k] + carry;
        long lo = prod & MASK32;
        long hi = prod >>> 32;

        long diff = uBuf[k] - lo;
        if (diff < 0) {
          diff += BASE;
          hi += 1;
        }

        U[j + i + k] = (int) diff;
        carry = hi;
      }
    }

    for (; i < n; i++) {
      long vi = VL[i];
      long ui = U[j + i] & MASK32;
      long prod = vi * qhat + carry;
      long lo = prod & MASK32;
      long hi = prod >>> 32;

      long diff = ui - lo;
      if (diff < 0) {
        diff += BASE;
        hi += 1;
      }

      U[j + i] = (int) diff;
      carry = hi;
    }

    long un = (U[j + n] & MASK32) - carry;
    U[j + n] = (int) un;

    return (un < 0) ? 1 : 0;
  }

  private static void addBack(final int[] U, final int j, final long[] VL, final int n, final int carryIn) {
    long carry = carryIn & MASK32;
    for (int i = 0; i < n; i++) {
      long sum = (U[i + j] & 0xFFFFFFFFL) + VL[i] + carry;
      U[i + j] = (int) sum;
      carry = sum >>> 32;
    }
    U[j + n] = (int) (U[j + n] + carry);
  }

  private static void addBackSIMD(final int[] U, final int j, final int[] V, final int n, final int carryIn) {
    final int lanes = ISPEC.length();
    long carry = carryIn & MASK32;
    int[] sumBuf = new int[lanes];

    int i = 0;
    for (; i + lanes <= n; i += lanes) {
      IntVector uVec = IntVector.fromArray(ISPEC, U, j + i);
      IntVector vVec = IntVector.fromArray(ISPEC, V, i);
      IntVector sVec = uVec.add(vVec);
      sVec.intoArray(sumBuf, 0);

      for (int k = 0; k < lanes; k++) {
        long sum = (sumBuf[k] & MASK32) + carry;
        if (sum >= BASE) {
          sum -= BASE;
          carry = 1;
        } else {
          carry = 0;
        }
          U[j + i + k] = (int) sum;
      }
    }

    for (; i < n; i++) {
      long sum = (U[j + i] & MASK32) + (V[i] & MASK32) + carry;
      if (sum >= BASE) {
        sum -= BASE;
        carry = 1;
      } else {
        carry = 0;
      }
      U[j + i] = (int) sum;
    }

    int idx = j + n;
    long sum = (U[idx] & MASK32) + carry;
    U[idx] = (int) (sum % BASE);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder("0x");
    for (byte b : toBytesBE()) {
      sb.append(String.format("%02x", b));
    }
    return sb.toString();
  }

  @Override
  public boolean equals(final Object obj) {
    if (this == obj) return true;
    if (!(obj instanceof UInt256)) return false;
    UInt256 other = (UInt256) obj;

    // Compare lengths after trimming leading zero limbs
    int cmp = UInt256.compare(this, other);
    return cmp == 0;
  }

  @Override
  public int hashCode() {
    int h = 1;
    for (int i = 0; i < length; i++) {
      h = 31 * h + limbs[i];
    }
    return h;
  }
}
