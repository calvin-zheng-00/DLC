import deeplabcut
import time

print (time.strftime("%H:%M:%S"))

deeplabcut.train_network('$PATH/config.yaml', displayiters=5000, saveiters=50000, maxiters=50000, keepdeconvweights=False)

deeplabcut.evaluate_network('$PATH/config.yaml', Shuffles=[1], plotting=True)

print (time.strftime("%H:%M:%S"))

