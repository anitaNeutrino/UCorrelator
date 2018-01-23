root -q doTMVA.C'(50,367,0)'
root -q doTMVA.C'(50,367,"simulated")'
root -q makeTMVATrees.C'("wais","wais",120,155,1,40)'
root -q makeTMVATrees.C'("a4all","anita4",50,367,40)'
root -q makeTMVATrees.C'("simulated","simulation",1,200,100)'
root -q fillTMVA.C
