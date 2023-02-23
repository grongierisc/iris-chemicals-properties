
CREATE MODEL pka_acidic PREDICTING (pka float) FROM Data.pka_acidic

TRAIN MODEL pka_acidic

VALIDATE MODEL pka_acidic

