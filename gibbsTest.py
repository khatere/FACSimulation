import gibbsEnergy

testTemp = 25.0

#Creating our gibbsEnergy Class
gibbsValues = gibbsEnergy.gibbs(testTemp)
y1 = gibbsValues.getDeltaGFe_Fe2()
y2 = gibbsValues.getDeltaGFe2Plus_Fe304()
y3 = gibbsValues.getDeltaGfH2_Hp()

print y1,y2,y3


gibbsValues.setTemperature(30.00)

y1 = gibbsValues.getDeltaGFe_Fe2()
y2 = gibbsValues.getDeltaGFe2Plus_Fe304()
y3 = gibbsValues.getDeltaGfH2_Hp()

print y1,y2,y3

