function Cell=DoG_WeightsOnly(Em)

[pL,pM,pS,cnAll,snAll,cpAll,spAll]=RandomizeConeInputs_Gauss2(Em);

Cell.CenterCones=[cnAll(1) cnAll(2) cnAll(3)];
Cell.SurroundCones=[snAll(1) snAll(2) snAll(3)];
Cell.CenterWeight=[cpAll(1) cpAll(2) cpAll(3)];
Cell.SurroundWeight=[spAll(1) spAll(2) spAll(3)];