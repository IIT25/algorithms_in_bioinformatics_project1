import pickle
import matplotlib.pyplot as plt
import pandas as pd

file = open("result.pkl",'rb')
results = pickle.load(file)
file.close()

results = pd.DataFrame(results)
results = results[results["size"] >10]

print(results)

plt.scatter(results["size"], results["time"])
plt.xlabel("n")
plt.ylabel("T(n)")
plt.show()

plt.scatter(results["size"], results["time"]/(results["size"]**2))
plt.xlabel("n")
plt.ylabel("T(n)/n^2")
plt.show()