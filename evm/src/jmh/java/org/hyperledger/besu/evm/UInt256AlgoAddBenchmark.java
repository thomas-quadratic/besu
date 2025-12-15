/*
 * Copyright contributors to Besu.
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

import java.math.BigInteger;
import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

@State(Scope.Thread)
@Warmup(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@Measurement(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@OutputTimeUnit(value = TimeUnit.NANOSECONDS)
@BenchmarkMode(Mode.AverageTime)
public class UInt256AlgoAddBenchmark {
  // Benchmarks for Addition

  protected static final int SAMPLE_SIZE = 30_000;

  protected byte[][] aPool;
  protected byte[][] bPool;
  protected int index;

  // Define available scenarios based on byte array sizes
  public enum Case {
    ADD_4_4(4, 4),        // 32-bit + 32-bit
    ADD_8_8(8, 8),        // 64-bit + 64-bit
    ADD_8_4(8, 4),        // 64-bit + 32-bit
    ADD_16_16(16, 16),    // 128-bit + 128-bit
    ADD_16_8(16, 8),      // 128-bit + 64-bit
    ADD_16_4(16, 4),      // 128-bit + 32-bit
    ADD_24_24(24, 24),    // 192-bit + 192-bit
    ADD_24_16(24, 16),    // 192-bit + 128-bit
    ADD_24_8(24, 8),      // 192-bit + 64-bit
    ADD_24_4(24, 4),      // 192-bit + 32-bit
    ADD_32_32(32, 32),    // 256-bit + 256-bit (full size)
    ADD_32_24(32, 24),    // 256-bit + 192-bit
    ADD_32_16(32, 16),    // 256-bit + 128-bit
    ADD_32_8(32, 8),      // 256-bit + 64-bit
    ADD_32_4(32, 4),      // 256-bit + 32-bit
    ADD_RAND(-1, -1);     // Random sizes

    final int aSize;
    final int bSize;

    Case(final int aSize, final int bSize) {
      this.aSize = aSize;
      this.bSize = bSize;
    }
  }

  @Param({
    "ADD_4_4",
    "ADD_8_8",
    "ADD_8_4",
    "ADD_16_16",
    "ADD_16_8",
    "ADD_16_4",
    "ADD_24_24",
    "ADD_24_16",
    "ADD_24_8",
    "ADD_24_4",
    "ADD_32_32",
    "ADD_32_24",
    "ADD_32_16",
    "ADD_32_8",
    "ADD_32_4",
    "ADD_RAND"
  })
  private String caseName;

  @Setup(Level.Iteration)
  public void setUp() {
    Case scenario = Case.valueOf(caseName);
    aPool = new byte[SAMPLE_SIZE][];
    bPool = new byte[SAMPLE_SIZE][];

    final ThreadLocalRandom random = ThreadLocalRandom.current();
    int aSize;
    int bSize;

    for (int i = 0; i < SAMPLE_SIZE; i++) {
      if (scenario.aSize < 0) {
        aSize = random.nextInt(1, 33);
      } else {
        aSize = scenario.aSize;
      }
      if (scenario.bSize < 0) {
        bSize = random.nextInt(1, 33);
      } else {
        bSize = scenario.bSize;
      }

      final byte[] a = new byte[aSize];
      final byte[] b = new byte[bSize];
      random.nextBytes(a);
      random.nextBytes(b);

      // Ensure positive numbers by clearing sign bit
      // if (a.length > 0) a[0] = (byte) (a[0] & 0x7F);
      // if (b.length > 0) b[0] = (byte) (b[0] & 0x7F);

      aPool[i] = a;
      bPool[i] = b;
    }
    index = 0;
  }

  // @Benchmark
  // public void addIntWidening(final Blackhole blackhole) {
  //   blackhole.consume(UInt256Algo.addIntWidening(aPool[index], bPool[index]));
  //   index = (index + 1) % SAMPLE_SIZE;
  // }

  // @Benchmark
  // public void addIntAndCarry(final Blackhole blackhole) {
  //   blackhole.consume(UInt256Algo.addIntAndCarry(aPool[index], bPool[index]));
  //   index = (index + 1) % SAMPLE_SIZE;
  // }

  // @Benchmark
  // public void addByteVarLen(final Blackhole blackhole) {
    // blackhole.consume(UInt256Algo.addByteVarLen(aPool[index], bPool[index]));
    // index = (index + 1) % SAMPLE_SIZE;
  // }

  // @Benchmark
  // public void addByteFixedLen(final Blackhole blackhole) {
    // blackhole.consume(UInt256Algo.addByteFixedLen(aPool[index], bPool[index]));
    // index = (index + 1) % SAMPLE_SIZE;
  // }

  // @Benchmark
  // public void addByteDualFixedLen(final Blackhole blackhole) {
    // blackhole.consume(UInt256Algo.addByteFixedLen(aPool[index], bPool[index]));
    // index = (index + 1) % SAMPLE_SIZE;
  // }

  @Benchmark
  public void addSIMDLong(final Blackhole blackhole) {
    blackhole.consume(UInt256Algo.addSIMDLong(aPool[index], bPool[index]));
    index = (index + 1) % SAMPLE_SIZE;
  }

  @Benchmark
  public void addSIMDInt(final Blackhole blackhole) {
    blackhole.consume(UInt256Algo.addSIMDInt(aPool[index], bPool[index]));
    index = (index + 1) % SAMPLE_SIZE;
  }
}
