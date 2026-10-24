# Cooperative Relay Transceiver — Part A

ICT 3128: transmitter design and implementation. This repository currently
contains Part A and its standalone functional validation.

## Run

Open this folder in MATLAB and run `partA_transmitter.m`. Keep all five `.m`
files together. Base MATLAB is sufficient; the code was run in MATLAB R2026b.
The demonstration regenerates the plots, allocation table, MAT data and
validation text in this folder.

## Files

- `transmitter.m`: reusable source function; call `tx = transmitter()`.
- `relay_decode_forward.m`: Part A relay detection, CRC and forwarding strategy.
- `qpsk_codec.m`: shared QPSK mapping and hard-decision demapping.
- `crc8_bits.m`: shared CRC-8 polynomial division.
- `partA_transmitter.m`: standalone demonstration and assertions.
- `partA_*.png`: payload, constellation, transmission timing and energy plots.
- `partA_power_allocation.csv`: configurable energy allocation table.
- `partA_validation.txt`: output of the verified MATLAB run.
- `partA_data.mat`: source and separate ideal/corrupted test structures.
- `Part_A_writeup.txt`: handwritten-report notes and integration contract.

## Defaults and interfaces

- Reproducibility: seed 16, twister; 100 frames.
- Frame: 16-bit ID (MSB first), 976 payload bits, 8-bit CRC; 500 QPSK symbols.
- CRC protects ID plus payload using x^8+x^2+x+1 (0x07), zero initialization,
  no reflection and no final XOR.
- QPSK: `[(1-2*b1) + j*(1-2*b2)]/sqrt(2)`.
- Energy budget: 500 normalized units per attempted frame, including both slots;
  default source/relay split 250/250. Symbol duration is normalized to one.
- Phase 1: source broadcasts. Phase 2: relay forwards on CRC pass, otherwise
  remains silent. Both time slots remain reserved on failure.
- Source and relay outputs already include their respective power scaling.
- Columns represent frames of symbol-level complex baseband data.
- The relay function takes coherent equalized symbols. The group must supply
  channel modeling and equalization before using actual received symbols.
- `sourceTwoPhase` and `relayTwoPhase` contain phase 1 followed by phase 2 in
  each frame column. Do not apply transmit-power scaling a second time.

Example configuration:

```matlab
tx = transmitter(struct('sourceEnergyFraction', 0.75));
relay = relay_decode_forward(tx.sourceTx, tx); % ideal input for a local test
```

## Validation and scope

The ideal local test recovered all 97,600 payload bits and complete frames
exactly, with 100/100 CRC passes. A separate one-bit corruption test rejected
frame 1 and forwarded the other 99 frames. This is a functional test, not a
measured wireless success rate. Per-frame budget checks passed for source
energy shares of 0.25, 0.50 and 0.75.

Relay forwarding is included in Part A as specified in the original assignment.
The group should reuse this implementation rather than duplicate it in Part C.
Channel/path-loss modeling, destination reception/combining and link-performance
analysis remain group integration work. Agree on the frame format and verify
exact noiseless recovery again with the group's receiver.

Direct and cooperative comparisons should use the same budget per attempted
payload frame. Energy per successfully delivered payload is an outcome metric.

The handwritten signed pages and consolidated group PDF are still required.
