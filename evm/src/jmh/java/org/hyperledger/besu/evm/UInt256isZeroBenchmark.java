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

import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.TimeUnit;

import org.apache.tuweni.bytes.Bytes;
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
@Warmup(iterations = 2, time = 1, timeUnit = TimeUnit.SECONDS)
@OutputTimeUnit(value = TimeUnit.NANOSECONDS)
@Measurement(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@BenchmarkMode(Mode.AverageTime)
public class UInt256isZeroBenchmark {
  protected static final int SAMPLE_SIZE = 30_000;
  protected UInt256[] pool;
  protected UInt256[] zPool;
  protected int index;

  @Setup(Level.Iteration)
  public void setUp() {
    pool = new UInt256[SAMPLE_SIZE];
    zPool = new UInt256[SAMPLE_SIZE];

    final ThreadLocalRandom random = ThreadLocalRandom.current();
    long[] limbs = new long[4];
    for (int i = 0; i < SAMPLE_SIZE; i++) {
      for (int j = 0; j < 4; j++) limbs[j] = random.nextLong();
      pool[i] = new UInt256(Arrays.copyOf(limbs, 4));
      for (int j = 0; j < 4; j++) limbs[j] = 0;
      zPool[i] = new UInt256(Arrays.copyOf(limbs, 4));
    }

    index = 0;
  }

  @Benchmark
  public void isZeroAgainstZero(final Blackhole blackhole) {
    blackhole.consume(zPool[index].isZero());
    index = (index + 1) % SAMPLE_SIZE;
  }

  @Benchmark
  public void isZeroArrayAgainstZero(final Blackhole blackhole) {
    blackhole.consume(zPool[index].isZeroArray());
    index = (index + 1) % SAMPLE_SIZE;
  }

  @Benchmark
  public void isZeroAgainstRandom(final Blackhole blackhole) {
    blackhole.consume(pool[index].isZero());
    index = (index + 1) % SAMPLE_SIZE;
  }

  @Benchmark
  public void isZeroArrayAgainstRandom(final Blackhole blackhole) {
    blackhole.consume(pool[index].isZeroArray());
    index = (index + 1) % SAMPLE_SIZE;
  }
}
