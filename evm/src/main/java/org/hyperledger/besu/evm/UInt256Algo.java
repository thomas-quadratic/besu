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
package org.hyperledger.besu.evm;

import java.util.Arrays;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.LongVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * 256-bits wide unsigned integer class.
 *
 * <p>This class is an optimised version of BigInteger for fixed width 256-bits integers.
 */
public final class UInt256Algo {
  // region Internals
  // --------------------------------------------------------------------------
  // UInt256 represents a big-endian 256-bits integer.
  // As opposed to Java int, operations are by default unsigned,
  // and signed version are interpreted in two-complements as usual.
  // Length is used to optimise algorithms, skipping leading zeroes.
  // Nonetheless, 256bits are always allocated and initialised to zeroes.

  /** Fixed sizes */
  public static final int BITSIZE = 256;
  public static final int BYTESIZE = 32;
  public static final int INTSIZE = 8;
  public static final int LONGSIZE = 4;

  // Fixed number of limbs or digits
  private static final int N_BITS_PER_LIMB = 32;
  // Mask for long values
  private static final long MASK_L = 0xFFFFFFFFL;

  // Arrays of zeros.
  private static final byte[] ZERO_BYTES = new byte[BYTESIZE];
  private static final byte[] ONE_BYTES =
      new byte[] {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
      };
  // For int, we accomodate up to a result of a multiplication
  private static final int[] ZERO_INTS = new int[INTSIZE + INTSIZE + 1];
  private static final long[] ZERO_LONGS = new long[LONGSIZE + LONGSIZE + 1];

  // SIMD Vector Species
  private static final VectorSpecies<Integer> SPECIES = IntVector.SPECIES_256;
  private static final VectorSpecies<Long> SPECIES_L = LongVector.SPECIES_256;
  // Vectors utils
  private static final LongVector ALLSET_VEC = LongVector.broadcast(SPECIES_L, -1L); // 0xFFFFFFFFFFFFFFFF

  // LookUp Table for reciprocal in long division
  private static final short[] LUT = new short[] {
        2045, 2037, 2029, 2021, 2013, 2005, 1998, 1990, 1983, 1975, 1968, 1960, 1953, 1946, 1938,
        1931, 1924, 1917, 1910, 1903, 1896, 1889, 1883, 1876, 1869, 1863, 1856, 1849, 1843, 1836,
        1830, 1824, 1817, 1811, 1805, 1799, 1792, 1786, 1780, 1774, 1768, 1762, 1756, 1750, 1745,
        1739, 1733, 1727, 1722, 1716, 1710, 1705, 1699, 1694, 1688, 1683, 1677, 1672, 1667, 1661,
        1656, 1651, 1646, 1641, 1636, 1630, 1625, 1620, 1615, 1610, 1605, 1600, 1596, 1591, 1586,
        1581, 1576, 1572, 1567, 1562, 1558, 1553, 1548, 1544, 1539, 1535, 1530, 1526, 1521, 1517,
        1513, 1508, 1504, 1500, 1495, 1491, 1487, 1483, 1478, 1474, 1470, 1466, 1462, 1458, 1454,
        1450, 1446, 1442, 1438, 1434, 1430, 1426, 1422, 1418, 1414, 1411, 1407, 1403, 1399, 1396,
        1392, 1388, 1384, 1381, 1377, 1374, 1370, 1366, 1363, 1359, 1356, 1352, 1349, 1345, 1342,
        1338, 1335, 1332, 1328, 1325, 1322, 1318, 1315, 1312, 1308, 1305, 1302, 1299, 1295, 1292,
        1289, 1286, 1283, 1280, 1276, 1273, 1270, 1267, 1264, 1261, 1258, 1255, 1252, 1249, 1246,
        1243, 1240, 1237, 1234, 1231, 1228, 1226, 1223, 1220, 1217, 1214, 1211, 1209, 1206, 1203,
        1200, 1197, 1195, 1192, 1189, 1187, 1184, 1181, 1179, 1176, 1173, 1171, 1168, 1165, 1163,
        1160, 1158, 1155, 1153, 1150, 1148, 1145, 1143, 1140, 1138, 1135, 1133, 1130, 1128, 1125,
        1123, 1121, 1118, 1116, 1113, 1111, 1109, 1106, 1104, 1102, 1099, 1097, 1095, 1092, 1090,
        1088, 1086, 1083, 1081, 1079, 1077, 1074, 1072, 1070, 1068, 1066, 1064, 1061, 1059, 1057,
        1055, 1053, 1051, 1049, 1047, 1044, 1042, 1040, 1038, 1036, 1034, 1032, 1030, 1028, 1026,
        1024,
  };

  // --------------------------------------------------------------------------
  // endregion

  // region Arithmetic Operations
  // --------------------------------------------------------------------------

  /**
   * Addition implem Int with Long widening
   *
   * @param a The left multi-precision integer
   * @param b The right multi-precision integer
   * @return The sum
   */
  public static byte[] addIntWidening(final byte[] a, final byte[] b) {
    if (isZero(a)) return b;
    if (isZero(b)) return a;
    int[] aLimbs = bytesToInts(a);
    int[] bLimbs = bytesToInts(b);
    int[] sum = new int[INTSIZE];
    addImplIntWidening(sum, aLimbs, bLimbs);
    return intsToBytes(sum, nLeadingZeroes(sum));
  }

  /**
   * Addition implem B
   *
   * @param a The left multi-precision integer
   * @param b The right multi-precision integer
   * @return The sum
   */
  public static byte[] addIntAndCarry(final byte[] a, final byte[] b) {
    if (isZero(a)) return b;
    if (isZero(b)) return a;
    int[] aLimbs = bytesToInts(a);
    int[] bLimbs = bytesToInts(b);
    int[] sum = new int[INTSIZE];
    addImplIntAndCarry(sum, aLimbs, bLimbs);
    return intsToBytes(sum, nLeadingZeroes(sum));
  }

  /**
   * Addition implem C
   *
   * @param a The left multi-precision integer
   * @param b The right multi-precision integer
   * @return The sum
   */
  public static byte[] addByteVarLen(final byte[] a, final byte[] b) {
    if (isZero(a)) return b;
    if (isZero(b)) return a;
    int len = Math.max(a.length, b.length);
    byte[] sum = new byte[len];
    byte carry = addImplByteVarLen(sum, a, b);
    byte[] result = sum;
    if (carry != 0) {
      result = new byte[sum.length + 1];
      result[0] = carry;
      System.arraycopy(sum, 0, result, 1, sum.length);
    }
    return result;
  }

  /**
   * Addition implementation D
   *
   * @param a The left multi-precision integer
   * @param b The right multi-precision integer
   * @return The sum
   */
  public static byte[] addByteFixedLen(final byte[] a, final byte[] b) {
    if (isZero(a)) return b;
    if (isZero(b)) return a;
    byte[] x;
    byte[] y;
    if (a.length < 32) {
      x = new byte[32];
      System.arraycopy(a, 0, x, 32 - a.length, a.length);
    } else {
      x = a;
    }
    if (b.length < 32) {
      y = new byte[32];
      System.arraycopy(b, 0, y, 32 - b.length, b.length);
    } else {
      y = b;
    }
    byte[] sum = new byte[BYTESIZE];
    addImplByteFixedLen(sum, x, y);
    return sum;
  }

  /**
   * Addition implementation E
   *
   * @param a The left multi-precision integer
   * @param b The right multi-precision integer
   * @return The sum
   */
  public static byte[] addByteDualFixedLen(final byte[] a, final byte[] b) {
    if (isZero(a)) return b;
    if (isZero(b)) return a;
    int len = Math.min(a.length, b.length);
    byte[] sum = new byte[BYTESIZE];
    byte carry = addImplByteDualFixedLen(sum, a, b, len);
    if (a.length > len) {
      carry = addImplByteDualFixedLen(sum, a, carry, len);
      len = a.length;
    } else if (b.length > len) {
      carry = addImplByteDualFixedLen(sum, b, carry, len);
      len = b.length;
    }
    if (carry != 0 && len < BYTESIZE) sum[BYTESIZE - len - 1] = carry;
    return sum;
  }

  public static byte[] addSIMDLong(final byte[] a, final byte[] b) {
    if (isZero(a)) return b;
    if (isZero(b)) return a;
    long[] aL = bytesToLongs(a);
    long[] bL = bytesToLongs(b);
    long[] sum = new long[LONGSIZE];
    addImplSIMDLong(sum, aL, bL);
    return longsToBytes(sum);
  }

  public static byte[] addSIMDInt(final byte[] a, final byte[] b) {
    if (isZero(a)) return b;
    if (isZero(b)) return a;
    int[] aI = bytesToInts(a);
    int[] bI = bytesToInts(b);
    int[] sum = new int[INTSIZE];
    addImplSIMDInt(sum, aI, bI);
    return intsToBytes(sum);
  }

  // ============================================================================
  // SIMD Addition Support
  // ============================================================================
  //
  // This section implements SIMD-accelerated 256-bit unsigned integer addition
  // using the JDK Vector API (jdk.incubator.vector). The implementation treats
  // 8 x 32-bit limbs as 4 x 64-bit lanes for hardware-accelerated parallel
  // addition.
  //
  // Hardware Requirements:
  // - x86-64 with AVX2 (256-bit SIMD) - widely available since 2013
  // - ARM with NEON/SVE (alternative vector instruction sets)
  //
  // Performance Benefits:
  // - 4 additions execute in parallel (SIMD lanes)
  // - Reduced pipeline stalls from data dependencies
  // - Efficient carry propagation via lookup table
  // - JIT compiler generates native SIMD instructions (e.g., vpaddd, vpaddq)
  //
  // Algorithm Overview:
  // 1. Load operands into 256-bit vectors (4 x 64-bit lanes)
  // 2. Parallel addition across all lanes simultaneously
  // 3. Detect per-lane carries using unsigned comparison (result < operand)
  // 4. Detect cascade conditions (lanes with all 1s propagate carries)
  // 5. Calculate cross-lane carry propagation (scalar operations)
  // 6. Apply cascaded carries via lookup table
  // 7. Return result with overflow flag
  //
  // Based on the reference implementation from .NET Runtime:
  // https://github.com/dotnet/runtime/blob/main/src/libraries/System.Runtime.Numerics/src/System/UInt256.cs
  // ============================================================================

  // Lookup table for SIMD carry propagation
  // Each entry represents carry values for 4 lanes based on cascade pattern (16 entries x 32 bytes)
  private static final byte[] BROADCAST_LOOKUP = {
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  };

  /**
   * SIMD-accelerated addition of two UInt256 values using JDK Vector API.
   *
   * <p>This implementation uses the JDK Vector API with 256-bit vectors (4 x 64-bit lanes) to
   * perform hardware-accelerated parallel addition with carry propagation across lanes.
   *
   * <p>Requires: --add-modules jdk.incubator.vector
   *
   * @param a First UInt256 operand
   * @param b Second UInt256 operand
   * @return Result UInt256 containing a + b, and whether overflow occurred
   */
  public static int addImplSIMDLong(final long[] result, final long[] a, final long[] b) {
    // Convert to 256-bit vectors with 4 x 64-bit lanes
    LongVector aVec = LongVector.fromArray(SPECIES_L, a, 0);
    LongVector bVec = LongVector.fromArray(SPECIES_L, b, 0);

    // Lanewise add, carry and cascade
    LongVector resultVec = aVec.add(bVec);
    VectorMask<Long> carryMask = resultVec.compare(VectorOperators.UNSIGNED_LT, aVec);
    VectorMask<Long> cascadeMask = resultVec.compare(VectorOperators.EQ, ALLSET_VEC);

    // Calculate cross-lane carries using scalar operations
    int carry = maskToInt(carryMask);
    int cascade = maskToInt(cascadeMask);
    carry = cascade + 2 * carry;
    cascade ^= carry;
    cascade &= 0x0f;

    int lookupOffset = cascade * 32;
    long[] cascadedCarries = new long[LONGSIZE];
    for (int i = 0; i < LONGSIZE; i++) {
      long value = 0;
      for (int j = 0; j < 8; j++) {
        value |= ((long) (BROADCAST_LOOKUP[lookupOffset + i * 8 + j] & 0xFF)) << (j * 8);
      }
      cascadedCarries[i] = value;
    }

    // Apply cascaded carries
    LongVector carryVector = LongVector.fromArray(SPECIES_L, cascadedCarries, 0);
    resultVec = resultVec.add(carryVector);
    resultVec.intoArray(result, 0);
    return ((carry & 0b1_0000) != 0) ? 1 : 0;
  }

  /**
   * Convert a VectorMask to an integer bitmask.
   *
   * <p>Each bit represents whether the corresponding lane's mask is true.
   *
   * @param mask The vector mask to convert
   * @return Integer bitmask where bit i is set if lane i's mask is true
   */
  private static int maskToInt(final VectorMask<Long> mask) {
    int result = 0;
    for (int i = 0; i < 4; i++) {
      if (mask.laneIsSet(i)) {
        result |= (1 << i);
      }
    }
    return result;
  }

  /**
   * SIMD-accelerated addition using IntVector (8 x 32-bit lanes).
   *
   * <p>This implementation uses IntVector for a more direct approach, performing parallel addition
   * across all 8 limbs, then handling carry propagation sequentially.
   *
   * <p>Simpler than the LongVector approach but requires sequential carry propagation.
   *
   * @param a First UInt256 operand
   * @param b Second UInt256 operand
   * @return Result UInt256 containing a + b, and whether overflow occurred
   */
  public static int addImplSIMDInt(final int[] result, final int[] a, final int[] b) {
    IntVector aVec = IntVector.fromArray(SPECIES, a, 0);
    IntVector bVec = IntVector.fromArray(SPECIES, b, 0);
    IntVector sumVec = aVec.add(bVec);
    VectorMask<Integer> carryMask = sumVec.compare(VectorOperators.UNSIGNED_LT, aVec);

    sumVec.intoArray(result, 0);
    long carry = 0;
    for (int i = 0; i < INTSIZE; i++) {
      if (carryMask.laneIsSet(i)) {
        carry = 1;
      }

      if (i < INTSIZE - 1 && carry != 0) {
        long next = (result[i + 1] & MASK_L) + carry;
        result[i + 1] = (int) next;
        carry = next >>> 32;
      }
    }
    return (int) carry;
  }

  /**
   * Multiplication
   *
   * @param a The left multi-precision integer
   * @param b The right multi-precision integer
   * @return The sum
   */
  public static byte[] mul(final byte[] a, final byte[] b) {
    if (isZero(a) || isZero(b)) return new byte[0];
    if (isOne(a)) return b;
    if (isOne(b)) return a;
    int[] aLimbs = bytesToInts(a, nLeadingZeroes(a));
    int[] bLimbs = bytesToInts(b, nLeadingZeroes(b));
    int[] prod = wrappingMul(aLimbs, 0, bLimbs, 0);
    return intsToBytes(prod, nLeadingZeroes(prod));
  }

  /**
   * Unsigned modulo reduction.
   *
   * @param a The dividend
   * @param b The modulus
   * @return The remainder of dividend modulo {@code modulus}.
   */
  public static byte[] mod(final byte[] a, final byte[] b) {
    if (isZero(a) || isZero(b)) return new byte[0];
    int[] aLimbs = bytesToInts(a, nLeadingZeroes(a));
    int[] bLimbs = bytesToInts(b, nLeadingZeroes(b));
    int[] remainder = knuthRemainder(aLimbs, bLimbs);
    return intsToBytes(remainder, nLeadingZeroes(remainder));
  }

  /**
   * Unsigned modulo reduction.
   *
   * @param a The dividend
   * @param b The modulus
   * @return The remainder of dividend modulo {@code modulus}.
   */
  public static byte[] modLongs(final byte[] a, final byte[] b) {
    if (isZero(a) || isZero(b)) return new byte[0];
    long[] aLimbs = bytesToLongs(a, nLeadingZeroes(a));
    long[] bLimbs = bytesToLongs(b, nLeadingZeroes(b));
    long[] remainder = knuthRemainder(aLimbs, bLimbs);
    return longsToBytes(remainder, nLeadingZeroes(remainder));
  }

  /**
   * Signed modulo reduction.
   *
   * <p>In signed modulo reduction, integers are interpretated as fixed 256 bits width two's
   * complement signed integers.
   *
   * @param a The dividend of the reduction
   * @param b The modulus of the reduction
   * @return The remainder modulo {@code modulus}.
   */
  public static byte[] signedMod(final byte[] a, final byte[] b) {
    if (isZero(a) || isZero(b)) return new byte[0];
    int[] x = bytesToInts(a);
    int[] y = bytesToInts(b);
    absInplace(x);
    absInplace(y);
    int[] r = knuthRemainder(x, y);
    if (isNeg(x)) {
      negate(r);
    }
    return intsToBytes(r);
  }

  /**
   * Modular addition.
   *
   * @param a The left dividend
   * @param b The right dividend
   * @param c The modulus of the reduction
   * @return a + b (mod c).
   */
  public static byte[] addMod(final byte[] a, final byte[] b, final byte[] c) {
    if (isZero(c) || isOne(c)) return new byte[0];
    if (isZero(a)) return mod(b, c);
    if (isZero(b)) return mod(a, c);
    int[] aLimbs = bytesToInts(a, nLeadingZeroes(a));
    int[] bLimbs = bytesToInts(b, nLeadingZeroes(b));
    int[] cLimbs = bytesToInts(c, nLeadingZeroes(c));
    int[] sum = new int[Math.max(aLimbs.length, bLimbs.length) + 1];
    addImplIntWidening(sum, aLimbs, bLimbs);
    int[] remainder = knuthRemainder(sum, cLimbs);
    return intsToBytes(remainder, nLeadingZeroes(remainder));
  }

  /**
   * Modular multiplication.
   *
   * @param a The left dividend
   * @param b The right dividend
   * @param c The modulus of the reduction
   * @return a * b (mod c).
   */
  public static byte[] mulMod(final byte[] a, final byte[] b, final byte[] c) {
    if (isZero(a) || isZero(b) || isZero(c) || isOne(c)) return new byte[0];
    if (isOne(a)) return mod(b, c);
    if (isOne(b)) return mod(a, c);
    int[] aLimbs = bytesToInts(a, nLeadingZeroes(a));
    int[] bLimbs = bytesToInts(b, nLeadingZeroes(b));
    int[] cLimbs = bytesToInts(c, nLeadingZeroes(c));
    int[] prod = addMul(aLimbs, 0, bLimbs, 0);
    int[] remainder = knuthRemainder(prod, cLimbs);
    return intsToBytes(remainder, nLeadingZeroes(remainder));
  }

  /**
   * Unsigned division
   *
   * @param a The dividend
   * @param b The divisor
   * @return a / b
   */
  public static byte[] div(final byte[] a, final byte[] b) {
    // TODO
    return a;
  }

  /**
   * Signed division
   *
   * @param a The dividend
   * @param b The divisor
   * @return a / b
   */
  public static byte[] sDiv(final byte[] a, final byte[] b) {
    // TODO
    return a;
  }

  // --------------------------------------------------------------------------
  // endregion

  // region Bitwise Operations
  // --------------------------------------------------------------------------

  /**
   * Bitwise AND operation
   *
   * @param a left byte array
   * @param b right byte array
   * @return bitwise AND of a and b.
   */
  public static byte[] and(final byte[] a, final byte[] b) {
    byte[] result = new byte[BYTESIZE];
    for (int i = 0; i < BYTESIZE; i++) {
      result[i] = (byte) (a[i] & b[i]);
    }
    return result;
  }

  /**
   * Bitwise XOR operation
   *
   * @param a left byte array
   * @param b right byte array
   * @return bitwise XOR of a and b.
   */
  public static byte[] xor(final byte[] a, final byte[] b) {
    byte[] result = new byte[BYTESIZE];
    for (int i = 0; i < BYTESIZE; i++) {
      result[i] = (byte) (a[i] ^ b[i]);
    }
    return result;
  }

  /**
   * Bitwise OR operation
   *
   * @param a left byte array
   * @param b right byte array
   * @return bitwise OR of a and b.
   */
  public static byte[] or(final byte[] a, final byte[] b) {
    byte[] result = new byte[BYTESIZE];
    for (int i = 0; i < BYTESIZE; i++) {
      result[i] = (byte) (a[i] | b[i]);
    }
    return result;
  }

  /**
   * Bitwise NOT operation
   *
   * @param a byte array
   * @return bitwise NOT of a.
   */
  public static byte[] not(final byte[] a) {
    byte[] result = new byte[BYTESIZE];
    for (int i = 0; i < BYTESIZE; i++) {
      result[i] = (byte) ~a[i];
    }
    return result;
  }

  // --------------------------------------------------------------------------
  // endregion

  // region : Helpers
  // --------------------------------------------------------------------------

  // private static byte[] trimLeadingZeroes(byte[] a) {
  //   int offset = nLeadingZeroes(a);
  //   byte[] result = new byte[a.length - offset];
  //   System.arraycopy(a, offset, result, 0, a.length - offset);
  //   return result;
  // }

  private static int nLeadingZeroes(final byte[] a) {
    return Arrays.mismatch(a, ZERO_BYTES); // Most significant byte index
  }

  private static int nLeadingZeroes(final int[] a) {
    return Arrays.mismatch(a, ZERO_INTS); // Most significant byte index
  }

  private static int nLeadingZeroes(final long[] a) {
    return Arrays.mismatch(a, ZERO_LONGS); // Most significant byte index
  }

  // private static int effectiveLength(final byte[] a) {
  //   return a.length - nLeadingZeroes(a);
  // }

  private static int effectiveLength(final int[] a) {
    return a.length - nLeadingZeroes(a);
  }

  /**
   * Instantiates a new UInt256 from byte array.
   *
   * @param bytes raw bytes in BigEndian order.
   * @return Big-endian UInt256 represented by the bytes.
   */
  private static int[] bytesToInts(final byte[] bytes) {
    return bytesToInts(bytes, 0);
  }

  /**
   * Instantiates a new UInt256 from byte array.
   *
   * @param bytes raw bytes in BigEndian order.
   * @return Big-endian UInt256 represented by the bytes.
   */
  private static int[] bytesToInts(final byte[] bytes, final int msb) {
    if (bytes.length == 0 || msb == -1 || msb >= bytes.length) return new int[INTSIZE];
    int[] limbs = new int[INTSIZE];
    int i = INTSIZE - 1; // Index in int array
    int b = bytes.length - 1; // Index in bytes array
    int limb;
    for (; b >= msb; i--) {
      int shift = 0;
      limb = 0;
      for (int j = 0; j < 4 && b >= msb; j++, b--, shift += 8) {
        limb |= ((bytes[b] & 0xFF) << shift);
      }
      limbs[i] = limb;
    }
    return limbs;
  }

  /**
   * Convert to BigEndian byte array.
   *
   * @return Big-endian ordered bytes for this UInt256 value.
   */
  private static byte[] intsToBytes(final int[] a) {
    return intsToBytes(a, 0);
  }

  /**
   * Convert to BigEndian byte array.
   *
   * @return Big-endian ordered bytes for this UInt256 value.
   */
  private static byte[] intsToBytes(final int[] a, final int offset) {
    byte[] result = new byte[4 * (a.length - offset)];
    for (int i = offset, j = 0; i < a.length; i++, j += 4) {
      int value = a[i];
      result[j] = (byte) (value >>> 24);
      result[j + 1] = (byte) (value >>> 16);
      result[j + 2] = (byte) (value >>> 8);
      result[j + 3] = (byte) value;
    }
    return result;
  }

  /**
   * Instantiates a new UInt256 from byte array.
   *
   * @param bytes raw bytes in BigEndian order.
   * @return Big-endian UInt256 represented by the bytes.
   */
  private static long[] bytesToLongs(final byte[] bytes) {
    return bytesToLongs(bytes, 0);
  }

  /**
   * Instantiates a new UInt256 from byte array.
   *
   * @param bytes raw bytes in BigEndian order.
   * @return Big-endian UInt256 represented by the bytes.
   */
  private static long[] bytesToLongs(final byte[] bytes, final int msb) {
    long[] limbs = new long[LONGSIZE];
    if (bytes.length == 0 || msb == -1 || msb >= bytes.length) return limbs;
    int i = LONGSIZE - 1; // Index in int array
    int b = bytes.length - 1; // Index in bytes array
    long limb;
    for (; b >= msb; i--) {
      int shift = 0;
      limb = 0;
      for (int j = 0; j < 8 && b >= msb; j++, b--, shift += 8) {
        limb |= ((bytes[b] & 0xFFL) << shift);
      }
      limbs[i] = limb;
    }
    return limbs;
  }

  /**
   * Convert to BigEndian byte array.
   *
   * @return Big-endian ordered bytes for this UInt256 value.
   */
  private static byte[] longsToBytes(final long[] a) {
    return longsToBytes(a, 0);
  }

  /**
   * Convert to BigEndian byte array.
   *
   * @return Big-endian ordered bytes for this UInt256 value.
   */
  private static byte[] longsToBytes(final long[] a, final int offset) {
    byte[] result = new byte[8 * (a.length - offset)];
    for (int i = offset, j = 0; i < a.length; i++, j += 8) {
      long value = a[i];
      result[j] = (byte) (value >>> 56);
      result[j + 1] = (byte) (value >>> 48);
      result[j + 2] = (byte) (value >>> 40);
      result[j + 3] = (byte) (value >>> 32);
      result[j + 4] = (byte) (value >>> 24);
      result[j + 5] = (byte) (value >>> 16);
      result[j + 6] = (byte) (value >>> 8);
      result[j + 7] = (byte) value;
    }
    return result;
  }

  // private static void intsIntoBytes(int[] in, int inFrom, byte[] out, int outFrom, int len) {
  //   for (int i = inFrom, j = outFrom; i < inFrom + len; i++, j += 4) {
  //     int value = in[i];
  //     out[j] = (byte) (value >>> 24);
  //     out[j + 1] = (byte) (value >>> 16);
  //     out[j + 2] = (byte) (value >>> 8);
  //     out[j + 3] = (byte) value;
  //   }
  // }

  // private static String bytesToHex(byte[] a) {
  //   StringBuilder sb = new StringBuilder("0x");
  //   for (byte b : a) {
  //     sb.append(String.format("%02x", b));
  //   }
  //   return sb.toString();
  // }

  // --------------------------------------------------------------------------
  // endregion

  // region : Comparisons
  // --------------------------------------------------------------------------

  /**
   * Is the value 0 ?
   *
   * @return true if this UInt256 value is 0.
   */
  private static boolean isZero(final byte[] bytes) {
    int msb = Arrays.mismatch(bytes, ZERO_BYTES);
    return (msb == -1 || msb >= bytes.length);
  }

  private static boolean isOne(final byte[] bytes) {
    int msb = Arrays.mismatch(bytes, 0, bytes.length, ONE_BYTES, BYTESIZE - bytes.length, BYTESIZE);
    return (msb == -1 || msb >= bytes.length);
  }

  /**
   * Compares two bytes array seen as big-endian unsigned integers.
   *
   * @param a left unsigned integer
   * @param b right unsigned integer
   * @return 0 if a == b, negative if a &lt; b and positive if a &gt; b.
   */
  public static int compare(final byte[] a, final byte[] b) {
    int ai = nLeadingZeroes(a);
    int bi = nLeadingZeroes(b);
    int cmp = Integer.compare(a.length - ai, b.length - bi);
    if (cmp != 0) return cmp;
    int i = Arrays.mismatch(a, ai, a.length, b, bi, b.length);
    return (i == -1 || i >= a.length - ai) ? 0 : Integer.compareUnsigned(a[i], b[i]);
  }

  // /**
  //  * Compares two UInt256.
  //  *
  //  * @param a left UInt256
  //  * @param b right UInt256
  //  * @return 0 if a == b, negative if a &lt; b and positive if a &gt; b.
  //  */
  // private static int compare(final int[] a, final int[] b) {
  //   int i = Arrays.mismatch(a, b);
  //   return (i == -1) ? 0 : Integer.compareUnsigned(a[i], b[i]);
  // }

  // --------------------------------------------------------------------------
  // endregion

  // region Support (private) Algorithms
  // --------------------------------------------------------------------------

  // Comparing two int subarrays as big-endian multi-precision integers.
  private static int compareLimbs(final int[] a, final int[] b) {
    if (a.length >= b.length) {
      int diffLen = a.length - b.length;
      int cmp = Arrays.mismatch(a, 0, diffLen, ZERO_INTS, 0, diffLen);
      if (cmp != -1) return 1;
      int i = Arrays.mismatch(a, diffLen, a.length, b, 0, b.length);
      return (i == -1) ? 0 : Integer.compareUnsigned(a[i + diffLen], b[i]);
    } else {
      int diffLen = b.length - a.length;
      int cmp = Arrays.mismatch(b, 0, diffLen, ZERO_INTS, 0, diffLen);
      if (cmp != -1) return -1;
      int i = Arrays.mismatch(a, 0, a.length, b, diffLen, b.length);
      return (i == -1) ? 0 : Integer.compareUnsigned(a[i], b[i + diffLen]);
    }
  }

  // Comparing two int subarrays as big-endian multi-precision integers.
  private static int compareLimbs(final long[] a, final long[] b) {
    if (a.length >= b.length) {
      int diffLen = a.length - b.length;
      int cmp = Arrays.mismatch(a, 0, diffLen, ZERO_LONGS, 0, diffLen);
      if (cmp != -1) return 1;
      int i = Arrays.mismatch(a, diffLen, a.length, b, 0, b.length);
      return (i == -1) ? 0 : Long.compareUnsigned(a[i + diffLen], b[i]);
    } else {
      int diffLen = b.length - a.length;
      int cmp = Arrays.mismatch(b, 0, diffLen, ZERO_LONGS, 0, diffLen);
      if (cmp != -1) return -1;
      int i = Arrays.mismatch(a, 0, a.length, b, diffLen, b.length);
      return (i == -1) ? 0 : Long.compareUnsigned(a[i], b[i + diffLen]);
    }
  }

  // Does two-complements represent a negative number: i.e. is leading bit set ?
  private static boolean isNeg(final int[] x) {
    return x[0] < 0;
  }

  // Negate in two-complements representation: bitwise NOT + 1
  // Inplace: modifies input x.
  private static void negate(final int[] x) {
    int carry = 1;
    for (int i = x.length - 1; i >= 0; i--) {
      x[i] = ~x[i] + carry;
      carry = (x[i] == 0 && carry == 1) ? 1 : 0;
    }
  }

  // Replaces x with its absolute value in two-complements representation
  private static void absInplace(final int[] x) {
    if (isNeg(x)) negate(x);
  }

  private static int shiftLeftInto(
      final int[] result, final int[] x, final int xOffset, final int shift) {
    // Unchecked: result should be initialised with zeroes
    // Unchecked: result length should be at least x.length + 1
    // Unchecked: 0 <= shift < N_BITS_PER_LIMB
    if (shift == 0) {
      int xLen = x.length - xOffset;
      int resultOffset = result.length - xLen;
      System.arraycopy(x, xOffset, result, resultOffset, xLen);
      return 0;
    }
    int carry = 0;
    int j = result.length - 1;
    for (int i = x.length - 1; i >= xOffset; i--, j--) {
      result[j] = (x[i] << shift) | carry;
      carry = x[i] >>> (N_BITS_PER_LIMB - shift);
    }
    return carry;
  }

  private static long shiftLeftInto(
      final long[] result, final long[] x, final int xOffset, final int shift) {
    // Unchecked: result should be initialised with zeroes
    // Unchecked: result length should be at least x.length + 1
    // Unchecked: 0 <= shift < N_BITS_PER_LIMB
    if (shift == 0) {
      int xLen = x.length - xOffset;
      int resultOffset = result.length - xLen;
      System.arraycopy(x, xOffset, result, resultOffset, xLen);
      return 0;
    }
    long carry = 0;
    int j = result.length - 1;
    for (int i = x.length - 1; i >= xOffset; i--, j--) {
      result[j] = (x[i] << shift) | carry;
      carry = x[i] >>> (64 - shift);
    }
    return carry;
  }

  private static int shiftRightInto(
      final int[] result, final int[] x, final int xOffset, final int shift) {
    // Unchecked: result length should be at least x.length
    // Unchecked: 0 <= shift < N_BITS_PER_LIMB
    if (shift == 0) {
      int xLen = x.length - xOffset;
      int resultOffset = result.length - xLen;
      System.arraycopy(x, xOffset, result, resultOffset, xLen);
      return 0;
    }
    int carry = 0;
    int j = result.length - x.length + xOffset;
    for (int i = xOffset; i < x.length; i++, j++) {
      result[j] = (x[i] >>> shift) | carry;
      carry = x[i] << (N_BITS_PER_LIMB - shift);
    }
    return carry;
  }

  private static long shiftRightInto(
      final long[] result, final long[] x, final int xOffset, final int shift) {
    // Unchecked: result length should be at least x.length
    // Unchecked: 0 <= shift < N_BITS_PER_LIMB
    if (shift == 0) {
      int xLen = x.length - xOffset;
      int resultOffset = result.length - xLen;
      System.arraycopy(x, xOffset, result, resultOffset, xLen);
      return 0;
    }
    long carry = 0;
    int j = result.length - x.length + xOffset;
    for (int i = xOffset; i < x.length; i++, j++) {
      result[j] = (x[i] >>> shift) | carry;
      carry = x[i] << (64 - shift);
    }
    return carry;
  }

  private static int addImplIntWidening(final int[] sum, final int[] a, final int[] b) {
    // Unchecked: result.length > INTSIZE
    // Unchecked: x.length == y.length == INTSIZE
    // Unchecked: INTSIZE == 8
    long carry = 0;
    int i = sum.length - 1;
    carry = adcL(sum, a[7], b[7], carry, i--);
    carry = adcL(sum, a[6], b[6], carry, i--);
    carry = adcL(sum, a[5], b[5], carry, i--);
    carry = adcL(sum, a[4], b[4], carry, i--);
    carry = adcL(sum, a[3], b[3], carry, i--);
    carry = adcL(sum, a[2], b[2], carry, i--);
    carry = adcL(sum, a[1], b[1], carry, i--);
    carry = adcL(sum, a[0], b[0], carry, i--);
    return (int) carry;
  }

  private static long adcL(
      final int[] sum, final int a, final int b, final long carry, final int index) {
    long aL = a & MASK_L;
    long bL = b & MASK_L;
    long s = aL + bL + carry;
    sum[index] = (int) s;
    return s >>> N_BITS_PER_LIMB;
  }

  private static int addImplIntAndCarry(final int[] sum, final int[] a, final int[] b) {
    // Unchecked both length 8
    int carry = 0;
    carry = adc(sum, a, b, carry, 7);
    carry = adc(sum, a, b, carry, 6);
    carry = adc(sum, a, b, carry, 5);
    carry = adc(sum, a, b, carry, 4);
    carry = adc(sum, a, b, carry, 3);
    carry = adc(sum, a, b, carry, 2);
    carry = adc(sum, a, b, carry, 1);
    carry = adc(sum, a, b, carry, 0);
    return carry;
  }

  private static int adc(
      final int[] sum, final int[] a, final int[] b, final int carryIn, final int i) {
    int s = a[i] + b[i];
    int carryOut = Integer.compareUnsigned(s, a[i]) < 0 ? 1 : 0;
    s += carryIn;
    carryOut |= (s == 0 && carryIn == 1) ? 1 : 0;
    sum[i] = s;
    return carryOut;
  }

  private static byte addImplByteVarLen(final byte[] result, final byte[] a, final byte[] b) {
    byte[] x;
    byte[] y;
    if (a.length < b.length) {
      x = b;
      y = a;
    } else {
      x = a;
      y = b;
    }
    int i = x.length - 1;
    int j = y.length - 1;
    int carry = 0;
    for (; j >= 0; i--, j--) {
      int xi = x[i] & 0xFF;
      int yi = y[j] & 0xFF;
      int sum = xi + yi + carry;
      result[i] = (byte) sum;
      carry = sum >>> 8;
    }
    for (; i >= 0 && carry != 0; i--) {
      int xi = x[i] & 0xFF;
      int sum = xi + carry;
      result[i] = (byte) sum;
      carry = sum >>> 8;
    }
    for (; i >= 0; i--) {
      result[i] = x[i];
    }
    return (byte) carry;
  }

  private static byte addImplByteFixedLen(final byte[] result, final byte[] a, final byte[] b) {
    int i = BYTESIZE - 1;
    int carry = 0;
    for (; i >= 0; i--) {
      int ai = a[i] & 0xFF;
      int bi = b[i] & 0xFF;
      int sum = ai + bi + carry;
      result[i] = (byte) sum;
      carry = sum >>> 8;
    }
    return (byte) carry;
  }

  private static byte addImplByteDualFixedLen(final byte[] result, final byte[] a, final byte[] b, final int len) {
    int i = a.length - 1;
    int j = b.length - 1;
    int k = result.length - 1;
    int carry = 0;
    for (; k >= result.length - len; i--, j--, k--) carry = adc(result, k, a[i], b[j], carry);
    return (byte) carry;
  }

  private static byte addImplByteDualFixedLen(final byte[] result, final byte[] a, final int carry, final int len) {
    int i = a.length - len - 1;
    int k = result.length - len - 1;
    int c = carry;
    for (; i >= 0; i--, k--) c = adc(result, k, a[i], c);
    return (byte) c;
  }

  private static int adc(final byte[] result, final int index, final byte a, final byte b, final int carry) {
    int ai = a & 0xFF;
    int bi = b & 0xFF;
    int sum = ai + bi + carry;
    result[index] = (byte) sum;
    return sum >>> 8;
  }

  private static int adc(final byte[] result, final int index, final byte a, final int carry) {
    int ai = a & 0xFF;
    int sum = ai + carry;
    result[index] = (byte) sum;
    return sum >>> 8;
  }

  private static int[] addMul(final int[] a, final int aOffset, final int[] b, final int bOffset) {
    // Shortest in outer loop, swap if needed
    int[] x;
    int[] y;
    int xOffset;
    int yOffset;
    if (a.length - aOffset < b.length - bOffset) {
      x = b;
      xOffset = bOffset;
      y = a;
      yOffset = aOffset;
    } else {
      x = a;
      xOffset = aOffset;
      y = b;
      yOffset = bOffset;
    }
    int[] lhs = new int[x.length + y.length - xOffset - yOffset];

    // Main algo
    int xLen = x.length - xOffset;
    for (int i = y.length - 1; i >= yOffset; i--) {
      long carry = 0;
      long yi = y[i] & MASK_L;

      int k = i + xLen - yOffset;
      for (int j = x.length - 1; j >= xOffset; j--, k--) {
        long prod = yi * (x[j] & MASK_L);
        long sum = (lhs[k] & MASK_L) + prod + carry;
        lhs[k] = (int) sum;
        carry = sum >>> N_BITS_PER_LIMB;
      }

      // propagate leftover carry
      while (carry != 0 && k >= 0) {
        long sum = (lhs[k] & MASK_L) + carry;
        lhs[k] = (int) sum;
        carry = sum >>> 32;
        k--;
      }
    }
    return lhs;
  }

  private static int[] wrappingMul(final int[] a, final int aOffset, final int[] b, final int bOffset) {
    // Shortest in outer loop, swap if needed
    int[] x;
    int[] y;
    int xOffset;
    int yOffset;
    if (a.length - aOffset < b.length - bOffset) {
      x = b;
      xOffset = bOffset;
      y = a;
      yOffset = aOffset;
    } else {
      x = a;
      xOffset = aOffset;
      y = b;
      yOffset = bOffset;
    }

    int xLen = x.length - xOffset;
    int yLen = y.length - yOffset;
    int maxLen = xLen + yLen + 1;
    int resLen = Math.min(INTSIZE, maxLen);
    int[] result = new int[resLen];

    for (int i = y.length - 1, m = resLen - 1; i >= yOffset && m >= 0; i--, m--) {
      long yi = y[i] & MASK_L;
      long carry = 0;
      int k = m;
      for (int j = x.length - 1; j >= xOffset && k >= 0; j--, k--) {
        long prod = yi * (x[j] & MASK_L);
        long sum = (result[k] & MASK_L) + prod + carry;
        result[k] = (int) sum;
        carry = sum >>> N_BITS_PER_LIMB;
      }

      // Propagate leftover carry, but only within INTSIZE bounds
      while (carry != 0 && k >= 0) {
        long sum = (result[k] & MASK_L) + carry;
        result[k] = (int) sum;
        carry = sum >>> N_BITS_PER_LIMB;
        k--;
      }
      // Any carry beyond position 0 is discarded (wrapping behavior)
    }
    return result;
  }
  // --------------------------------------------------------------------------
  // endregion

  // region Division
  // ---------------------------------------------------------------------------
  private static int[] knuthRemainder(final int[] dividend, final int[] modulus) {
    // Unchecked: modulus is non Zero and non One.
    int[] result = new int[INTSIZE];
    int modLen = effectiveLength(modulus);
    int divLen = effectiveLength(dividend);

    // Shortcut: if dividend < modulus or dividend == modulus
    int cmp = compareLimbs(dividend, modulus);
    if (cmp < 0) {
      System.arraycopy(dividend, dividend.length - divLen, result, INTSIZE - divLen, divLen);
      return result;
    } else if (cmp == 0) {
      return result;
    }

    // Shortcut: if modulus has a single limb
    if (modLen == 1) {
      if (divLen == 1) {
        result[INTSIZE - 1] =
            Integer.remainderUnsigned(dividend[dividend.length - 1], modulus[modulus.length - 1]);
        return result;
      }
      long d = modulus[modulus.length - 1] & MASK_L;
      long rem = 0;
      // Process from most significant limb downwards
      for (int i = dividend.length - divLen; i < dividend.length; i++) {
        long cur = (rem << 32) | (dividend[i] & MASK_L);
        rem = Long.remainderUnsigned(cur, d);
      }
      result[INTSIZE - 1] = (int) rem;
      result[INTSIZE - 2] = (int) (rem >>> 32);
      return result;
    }

    int shift = Integer.numberOfLeadingZeros(modulus[modulus.length - modLen]);
    // Normalize
    int[] vLimbs = new int[modLen];
    shiftLeftInto(vLimbs, modulus, modulus.length - modLen, shift);
    int[] uLimbs = new int[divLen + 1];
    uLimbs[0] = shiftLeftInto(uLimbs, dividend, dividend.length - divLen, shift);
    int diffLen = divLen - modLen + 1;

    long[] vLimbsAsLong = new long[modLen];
    for (int i = 0; i < modLen; i++) {
      vLimbsAsLong[i] = vLimbs[i] & MASK_L;
    }

    // Main division loop
    long vn1 = vLimbsAsLong[0];
    long vn2 = vLimbsAsLong[1];
    for (int j = 1; j < diffLen + 1; j++) {
      long ujn = (uLimbs[j - 1] & MASK_L);
      long ujn1 = (uLimbs[j] & MASK_L);
      long ujn2 = (uLimbs[j + 1] & MASK_L);

      long dividendPart = (ujn << N_BITS_PER_LIMB) | ujn1;
      // Check that no need for Unsigned version of divrem.
      long qhat = Long.divideUnsigned(dividendPart, vn1);
      long rhat = Long.remainderUnsigned(dividendPart, vn1);

      while (qhat == 0x1_0000_0000L
          || Long.compareUnsigned(qhat * vn2, (rhat << N_BITS_PER_LIMB) | ujn2) > 0) {
        qhat--;
        rhat += vn1;
        if (rhat >= 0x1_0000_0000L) break;
      }

      // Multiply-subtract qhat*v from u slice
      long borrow = 0;
      for (int i = modLen - 1; i >= 0; i--) {
        long prod = vLimbsAsLong[i] * qhat;
        long sub = (uLimbs[i + j] & MASK_L) - (prod & MASK_L) - borrow;
        uLimbs[i + j] = (int) sub;
        borrow = (prod >>> N_BITS_PER_LIMB) - (sub >> N_BITS_PER_LIMB);
      }
      long sub = (uLimbs[j - 1] & MASK_L) - borrow;
      uLimbs[j - 1] = (int) sub;

      if (sub < 0) {
        // Add back
        long carry = 0;
        for (int i = modLen - 1; i >= 0; i--) {
          long sum = (uLimbs[i + j] & MASK_L) + vLimbsAsLong[i] + carry;
          uLimbs[i + j] = (int) sum;
          carry = sum >>> N_BITS_PER_LIMB;
        }
        uLimbs[j - 1] = (int) (uLimbs[j - 1] + carry);
      }
    }
    // Unnormalize remainder
    shiftRightInto(result, uLimbs, diffLen, shift);
    return result;
  }

  private static long reciprocal(final long x) {
    // Unchecked: x >= (1 << 63)
    long x0 = x & 1L;
    int x9 = (int) (x >>> 55);
    long x40 = 1 + (x >>> 24);
    long x63 = (x + 1) >>> 1;
    long v0 = LUT[x9 - 256] & 0xFFFFL;
    long v1 = (v0 << 11) - ((v0 * v0 * x40) >>> 40) - 1;
    long v2 = (v1 << 13) + ((v1 * ((1L << 60) - v1 * x40)) >>> 47);
    long e = ((v2 >>> 1) & (-x0)) - v2 * x63;
    long s = Math.unsignedMultiplyHigh(v2, e);
    long v3 = (s >>> 1) + (v2 << 31);
    long t0 = v3 * x;
    long t1 = Math.unsignedMultiplyHigh(v3, x);
    t0 += x;
    if (Long.compareUnsigned(t0, x) < 0) t1++;
    t1 += x;
    long v4 = v3 - t1;
    return v4;
  }

  private static long reciprocal(final long x0, final long x1) {
    // Unchecked: <x1, x0> >= (1 << 127)
    long v = reciprocal(x1);
    long p = x1 * v + x0;
    if (Long.compareUnsigned(p, x0) < 0) {
      v--;
      if (Long.compareUnsigned(p, x1) >= 0) {
        v--;
        p -= x1;
      }
      p -= x1;
    }
    long t0 = v * x0;
    long t1 = Math.unsignedMultiplyHigh(v, x0);
    p += t1;
    if (Long.compareUnsigned(p, t1) < 0) {
      v--;
      int cmp = Long.compareUnsigned(p, x1);
      if ((cmp > 0) || ((cmp == 0) && (Long.compareUnsigned(t0, x0) >= 0))) v--;
    }
    return v;
  }

  private static long div2by1(final long[] xs, final int xOffset, final long y, final long yInv) {
    // Performs on x a single step of x / y
    long x0 = xs[xOffset + 1];
    long x1 = xs[xOffset];

    // wrapping umul x1 * yInv
    long q0 = x1 * yInv;
    long q1 = Math.unsignedMultiplyHigh(x1, yInv);

    // wrapping uadd <q1, q0> + <x1, x0> + <1, 0>
    long sum = q0 + x0;
    long carry = ((q0 & x0) | ((q0 | x0) & ~sum)) >>> 63;
    q0 = sum;
    q1 += x1 + carry + 1;

    xs[xOffset + 1] -= q1 * y;

    if (Long.compareUnsigned(xs[xOffset + 1], q0) > 0) {
      q1 -= 1;
      xs[xOffset + 1] += y;
    }

    if (Long.compareUnsigned(xs[xOffset + 1], y) >= 0) {
      q1 += 1;
      xs[xOffset + 1] -= y;
    }
    xs[xOffset] = 0;
    return q1;
  }

  private static long div3by2(
      final long[] xs, final int xOffset, final long y0, final long y1, final long yInv) {
    // <x2, x1, x0> divided by <y1, y0>.
    // Works inplace on x, modifying it to hold the remainder.
    // Returns the quotient q.
    // Requires <x2, x1> < <y1, y0> otherwise quotient overflows.
    long overflow; // carry or borrow
    long res;  // sum or diff
    long x0 = xs[xOffset + 2];
    long x1 = xs[xOffset + 1];
    long x2 = xs[xOffset];

    // <q1, q0> = x2 * yInv + <x2, x1>
    long q0 = x2 * yInv;
    long q1 = Math.unsignedMultiplyHigh(x2, yInv);
    res = q0 + x1;
    overflow = ((q0 & x1) | ((q0 | x1) & ~res)) >>> 63;
    q0 = res;
    q1 = q1 + x2 + overflow;

    // r1 <- x1 - q1 * y1 mod B
    xs[xOffset + 1] -= q1 * y1;

    // wrapping sub <r1, x0> − q1*y0 − <y1, y0>
    long t0 = q1 * y0;
    long t1 = Math.unsignedMultiplyHigh(q1, y0);

    res = x0 - t0;
    overflow = ((~x0 & t0) | ((~x0 | t0) & res)) >>> 63;
    xs[xOffset + 2] = res;
    xs[xOffset + 1] -= (t1 + overflow);

    res = xs[xOffset + 2] - y0;
    overflow = ((~xs[xOffset + 2] & y0) | ((~xs[xOffset + 2] | y0) & res)) >>> 63;
    xs[xOffset + 2] = res;
    xs[xOffset + 1] -= (y1 + overflow);

    // Adjustments
    q1 += 1;
    if (Long.compareUnsigned(xs[xOffset + 1], q0) >= 0) {
      q1 -= 1;
      res = xs[xOffset + 2] + y0;
      overflow = ((xs[xOffset + 2] & y0) | ((xs[xOffset + 2] | y0) & ~res)) >>> 63;
      xs[xOffset + 2] = res;
      xs[xOffset + 1] += y1 + overflow;
    }

    int cmp = Long.compareUnsigned(xs[xOffset + 1], y1);
    if ((cmp > 0) || ((cmp == 0) && (Long.compareUnsigned(xs[xOffset + 2], y0) >= 0))) {
      q1 += 1;
      res = xs[xOffset + 2] - y0;
      overflow = ((~xs[xOffset + 2] & y0) | ((~xs[xOffset + 2] | y0) & res)) >>> 63;
      xs[xOffset + 2] = res;
      xs[xOffset + 1] -= (y1 + overflow);
    }
    xs[xOffset] = 0;
    return q1;
  }

  private static void modNby1(final long[] xs, final int xOffset, final long y, final long yInv) {
    // Unchecked: y > 1 << 63
    // if xs >= y, overflows, so must substract y first
    // Can happen only with most significant digit ?
    for (int i = xOffset; i < xs.length - 1; i++) {
      if (Long.compareUnsigned(xs[i], y) >= 0) xs[i] -= y;
      div2by1(xs, i, y, yInv);
    }
  }

  private static void modNby2(
      final long[] xs, final int xOffset, final long y0, final long y1, final long yInv) {
    // if <x+1, x+0> >= <y1, y0>, overflows, so must substract <y1, y0> first
    for (int i = xOffset; i < xs.length - 2; i++) {
      long x1 = xs[i];
      long x0 = xs[i + 1];
      int cmp = Long.compareUnsigned(x1, y1);
      if ((cmp > 0) || ((cmp == 0) && (Long.compareUnsigned(x0, y0) >= 0))) {
        long borrow = (Long.compareUnsigned(y0, x0) > 0) ? 1 : 0;
        xs[i + 1] -= y0;
        xs[i] -= (y1 + borrow);
      }
      long q = div3by2(xs, i, y0, y1, yInv);
    }
  }

  private static void modNbyM(
      final long[] xs, final int xOffset, final long[] ys, final int yOffset, final long yInv) {
    long res;
    long borrow = 0;
    long y1 = ys[yOffset];
    long y0 = ys[yOffset + 1];
    int yLen = ys.length - yOffset;
    int xLen = xs.length - xOffset;
    for (int j = 0; j < xLen - yLen; j++) {
      // Unlikely overflow 3x2 case: x[..2] == y
      if (((y1 ^ xs[j]) | (y0 ^ xs[j + 1])) == 0) {
        // q = <1, 0> in this case, so multiply is trivial.
        // Still need to substract (and add back if needed).
        int i = j + yLen;
        for (int k = ys.length - 1; k >= yOffset; i--, k--) {
          res = xs[i] - ys[k] - borrow;
          borrow = ((~xs[i] & ys[k]) | ((~xs[i] | ys[k]) & res)) >>> 63;
          xs[i] = res;
        }
        // xs[i] -= borrow;
        // borrow = (xs[i] >>> 63) & 1L;
      } else {
        long q = div3by2(xs, j, y0, y1, yInv);
        long res0;
        long prodHigh = 0;

        // Multiply-subtract: already have highest 2 limbs
        int i = j + yLen;
        for (int k = ys.length - 1; k >= yOffset + 2; k--, i--) {
          long p0 = ys[k] * q;
          long p1 = Math.unsignedMultiplyHigh(ys[k], q);
          res0 = xs[i] - p0 - borrow;
          p1 += ((~xs[i] & p0) | ((~xs[i] | p0) & res0)) >>> 63;  // p1 < b - 1 so can add 0 or 1
          res = res0 - prodHigh;
          borrow = ((~res0 & prodHigh) | ((~res0 | prodHigh) & res)) >>> 63;
          xs[i] = res;
          prodHigh = p1;
        }
        // Propagate overflows (borrows)
        res = xs[i] - prodHigh - borrow;
        borrow = ((~xs[i] & prodHigh) | ((~xs[i] | prodHigh) & res)) >>> 63;
        xs[i--] = res;

        res = xs[i] - borrow;
        borrow = (((borrow ^ 1L) | xs[i]) == 0) ? 1 : 0;
        xs[i--] = res;

        xs[i] -= borrow;
        borrow = (xs[i] >>> 63) & 1;
      }

      if (borrow != 0) { // unlikely
        // Add back
        long carry = 0;
        int i = j + yLen;
        for (int k = ys.length - 1; k >= yOffset; k--, i--) {
          res = xs[i] + ys[k] + carry;
          carry = ((xs[i] & ys[k]) | ((xs[i] | ys[k]) & ~res)) >>> 63;
          xs[i] = res;
        }
        for (; (carry != 0) && (i >= j); i--) {
          xs[i] += carry;
          carry = (xs[i] == 0L) ? 1L : 0L;
        }
      }
      borrow = 0;
    }
  }

  private static long[] knuthRemainder(final long[] dividend, final long[] modulus) {
    long[] result = new long[LONGSIZE];
    int divLen = dividend.length - nLeadingZeroes(dividend);
    int modLen = modulus.length - nLeadingZeroes(modulus);


    // Shortcuts
    if (modLen == 0) return result;
    if (divLen == 1) {
      if (modLen == 1) result[result.length - 1] = Long.remainderUnsigned(dividend[dividend.length - 1], modulus[dividend.length - 1]);
      else System.arraycopy(dividend, dividend.length - divLen, result, LONGSIZE - divLen, divLen);
      return result;
    }
    // Shortcut: if dividend < modulus or dividend == modulus
    int cmp = compareLimbs(dividend, modulus);
    if (cmp < 0) {
      System.arraycopy(dividend, dividend.length - divLen, result, LONGSIZE - divLen, divLen);
      return result;
    } else if (cmp == 0) {
      return result;
    }

    // Perform Division
    // -- Normalise
    int shift = Long.numberOfLeadingZeros(modulus[modulus.length - modLen]);
    long[] vLimbs = new long[modLen];
    shiftLeftInto(vLimbs, modulus, modulus.length - modLen, shift);
    long[] uLimbs = new long[divLen + 1];
    uLimbs[0] = shiftLeftInto(uLimbs, dividend, dividend.length - divLen, shift);
    int diffLen = divLen - modLen + 1;

    // -- Divide
    if (modLen == 1) {
      long inv = reciprocal(vLimbs[0]);
      modNby1(uLimbs, 0, vLimbs[0], inv);
    } else if (modLen == 2) {
      long inv = reciprocal(vLimbs[1], vLimbs[0]);
      modNby2(uLimbs, 0, vLimbs[1], vLimbs[0], inv);
    } else {
      long inv = reciprocal(vLimbs[modLen - 2], vLimbs[modLen - 1]);
      modNbyM(uLimbs, 0, vLimbs, modLen, inv);
    }
    // -- Unnormalize
    shiftRightInto(result, uLimbs, diffLen, shift);
    return result;
  }
  // --------------------------------------------------------------------------
  // endregion
}
