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
package org.hyperledger.besu.datatypes;

import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.TimeUnit;

import org.apache.tuweni.bytes.Bytes;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

@State(Scope.Thread)
@Warmup(iterations = 2, time = 1, timeUnit = TimeUnit.SECONDS)
@OutputTimeUnit(value = TimeUnit.NANOSECONDS)
@Measurement(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@BenchmarkMode(Mode.AverageTime)
public class UInt256Benchmark {

  protected static final int SAMPLE_SIZE = 30_000;

  protected byte[][] bytesPool;
  protected UInt256[] uint256Pool;
  protected int[] shiftPool;
  protected int index;

  @Setup()
  public void setUp() {
    final ThreadLocalRandom random = ThreadLocalRandom.current();
    bytesPool = new byte[SAMPLE_SIZE][];
    uint256Pool = new UInt256[SAMPLE_SIZE];
    shiftPool = new int[SAMPLE_SIZE];
    for (int i = 0; i < SAMPLE_SIZE; i++) {
      final int size = 1 + random.nextInt(32); // [1, 32]
      final byte[] bytes = new byte[size];
      random.nextBytes(bytes);
      bytesPool[i] = bytes;
      uint256Pool[i] = UInt256.fromBytesBE(bytes);
      shiftPool[i] = 1 + random.nextInt(Math.max(1, size - 1)); // [1, size-1]
    }
    index = 0;
  }

  @Benchmark
  public void baseline(final Blackhole blackhole) {
    final int i = index;
    index = (index + 1) % SAMPLE_SIZE;
    blackhole.consume(i);
  }

  @Benchmark
  public void fromBytesBE(final Blackhole blackhole) {
    final int idx = index;
    index = (index + 1) % SAMPLE_SIZE;
    blackhole.consume(UInt256.fromBytesBE(bytesPool[idx]));
  }

  @Benchmark
  public void shiftLeft(final Blackhole blackhole) {
    final int i = index;
    index = (index + 1) % SAMPLE_SIZE;
    blackhole.consume(uint256Pool[i].shiftLeft(shiftPool[i]));

  }

  @Benchmark
  public void shiftRight(final Blackhole blackhole) {
    final int i = index;
    index = (index + 1) % SAMPLE_SIZE;
    blackhole.consume(uint256Pool[i].shiftRight(shiftPool[i]));
  }

  @Benchmark
  public void shiftTooMuch(final Blackhole blackhole) {
    final int i = index;
    index = (index + 1) % SAMPLE_SIZE;
    blackhole.consume(uint256Pool[i].shiftLeft(35));
  }

  @Benchmark
  public void shiftNoBitShift(final Blackhole blackhole) {
    final int i = index;
    index = (index + 1) % SAMPLE_SIZE;
    blackhole.consume(uint256Pool[i].shiftLeft(4));
  }
}
