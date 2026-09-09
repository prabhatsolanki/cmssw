#!/usr/bin/env python3
"""
export_vertex_slot_for_alpaka.py - wrap a vertex-slot TorchScript model for the Alpaka (PyTorchAlpaka) backend

The training export (VertexSlotModelInference) returns A[N,K], z_hat[K], t_hat[K], p[K], pi[N,3].
PyTorchAlpaka maps every output to a SoA column with one row per track, so z_hat, t_hat and p are
replicated to [N,K]. K must be equal to vertexgnn::kNumSlots (DataFormats/VertexGNNReco/interface/VertexGNNSoA.h)
and the number of input features to the TrackFeaturesLayout columns of the same header.

Usage:  python3 export_vertex_slot_for_alpaka.py <model.pt> <out_alpaka.pt> --num-slots 220 --num-features 13
"""
import argparse
from typing import Tuple

import torch
import torch.nn as nn


class AlpakaWrapper(nn.Module):
    def __init__(self, inner: torch.jit.ScriptModule):
        super().__init__()
        self.inner = inner

    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        A, z_hat, t_hat, p, pi = self.inner(x)
        N = A.size(0)
        return A, z_hat.unsqueeze(0).expand(N, -1), t_hat.unsqueeze(0).expand(N, -1), p.unsqueeze(0).expand(N, -1), pi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model"); ap.add_argument("output")
    ap.add_argument("--num-slots", type=int, required=True, help="kNumSlots of the CMSSW build")
    ap.add_argument("--num-features", type=int, required=True, help="number of TrackFeaturesLayout columns")
    args = ap.parse_args()

    original = torch.jit.load(args.model, map_location="cpu").eval()
    n = 100
    x = torch.randn(1, n, args.num_features)
    x[0, :, 1] = x[0, :, 1].abs() + 0.01       # dz > 0
    x[0, :, 9:12] = x[0, :, 9:12].abs() + 0.01  # sigma_t > 0
    x[0, :, 12] = (torch.rand(n) > 0.3).float()  # has_time
    with torch.no_grad():
        A, z_hat, t_hat, p, pi = original(x)
    K = z_hat.shape[0]
    print(f"model: A{tuple(A.shape)} z_hat{tuple(z_hat.shape)} t_hat{tuple(t_hat.shape)} p{tuple(p.shape)} pi{tuple(pi.shape)} -> K={K}")
    if K != args.num_slots:
        raise SystemExit(f"ERROR: model has K={K} slots but the CMSSW build expects kNumSlots={args.num_slots}")

    wrapped = torch.jit.script(AlpakaWrapper(original).eval())
    wrapped.save(args.output)
    loaded = torch.jit.load(args.output, map_location="cpu")
    with torch.no_grad():
        Aw, zw, tw, pw, piw = loaded(x)
        A2, z2, t2, p2, pi2 = loaded(x.squeeze(0))  # CMSSW passes [N, features]
    worst = max((A - Aw).abs().max().item(), (z_hat - zw[0]).abs().max().item(), (t_hat - tw[0]).abs().max().item(),
                (p - pw[0]).abs().max().item(), (pi - piw).abs().max().item(), (Aw - A2).abs().max().item(),
                (zw - z2).abs().max().item(), (pw - p2).abs().max().item())
    print(f"wrapped outputs: A{tuple(Aw.shape)} z_hat{tuple(zw.shape)} t_hat{tuple(tw.shape)} p{tuple(pw.shape)} pi{tuple(piw.shape)}; max |diff| vs original = {worst:.2e}")
    if worst > 1e-6:
        raise SystemExit("ERROR: wrapped model does not reproduce the original")
    print(f"saved {args.output}")


if __name__ == "__main__":
    main()
