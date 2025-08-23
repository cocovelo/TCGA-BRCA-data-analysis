import matplotlib.pyplot as plt
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
import pandas as pd
from sklearn.preprocessing import StandardScaler
import numpy as np

df_expression = pd.read_csv("C:/Users/colin/PycharmProjects/PythonProject/tcga-brca-project/tcga-brca-vst-normalized-counts.csv", index_col=0)
meta_df = pd.read_csv("C:/Users/colin/PycharmProjects/PythonProject/tcga-brca-project/tcga_brca_sample_metadata.csv", index_col=0)

X = df_expression.T
common_samples = list(set(X.index) & set(meta_df.index))
X = X.loc[common_samples].sort_index()
meta_df = meta_df.loc[common_samples].sort_index()

scaler = StandardScaler()
X_scaled = pd.DataFrame(scaler.fit_transform(X), columns=X.columns, index=X.index)

target_variable = 'paper_BRCA_Subtype_PAM50'
if target_variable not in meta_df.columns:
    if 'vital_status' in meta_df.columns:
        target_variable = 'vital_status'
        print(f"Using '{target_variable}' as target variable instead.")
    else:
        raise ValueError("No suitable target variable found for kNN classification.")

y_labels = meta_df[target_variable]
valid_indices = y_labels.dropna().index
X_knn = X_scaled.loc[valid_indices]
y_knn = y_labels.loc[valid_indices]

label_encoder = LabelEncoder()
y_encoded = label_encoder.fit_transform(y_knn)

X_train, X_test, y_train, y_test = train_test_split(
    X_knn, y_encoded, test_size=0.2, random_state=42, stratify=y_encoded
)

k_values = range(1, 31)
accuracies = []

for k in k_values:
    knn = KNeighborsClassifier(n_neighbors=k)
    knn.fit(X_train, y_train)
    y_pred = knn.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    accuracies.append(accuracy)

plt.figure(figsize=(10, 6))
plt.plot(k_values, accuracies, marker='o', linestyle='-', color='skyblue')
plt.title('kNN Classifier Accuracy vs. Number of Neighbors (k)')
plt.xlabel('Number of Neighbors (k)')
plt.ylabel('Accuracy on Test Set')
plt.xticks(np.arange(min(k_values), max(k_values)+1, 2))
plt.grid(True, linestyle='--', alpha=0.7)
plt.axvline(x=k_values[np.argmax(accuracies)], color='red', linestyle='--', label=f'Optimal k: {k_values[np.argmax(accuracies)]}')
plt.legend()
plt.tight_layout()
plt.show()

print(f"\nAccuracies for k values {min(k_values)}-{max(k_values)}: {accuracies}")
optimal_k_index = np.argmax(accuracies)
optimal_k = k_values[optimal_k_index]
max_accuracy = accuracies[optimal_k_index]
print(f"The optimal k value is {optimal_k} with an accuracy of {max_accuracy:.4f}.")