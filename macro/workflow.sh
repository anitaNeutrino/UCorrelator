root -q doTMVA.C'(40,367,0)'
root -q doTMVA.C'(40,367,"simulated")'
root -q makeTMVATrees.C'("wais","wais",120,155,1,40)'
root -q makeTMVATrees.C'("a4all","anita4",41,367,40)'
root -q makeTMVATrees.C'("simulated","simulation",1,500,40)'
root -q fillTMVA.C
root plotFisher.C
root -q makeSourceMap.C'("wais")'
root -q makeSourceMap.C'("anita4")'
root -q makeSourceMap.C'("simulation")'
