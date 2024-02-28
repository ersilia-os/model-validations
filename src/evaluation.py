import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix, roc_curve, roc_auc_score

def calculate_metrics_and_display(predictions, true_labels, model_name):
    predictions_positive_class = predictions['probability']
    predictions_positive_class = np.nan_to_num(predictions_positive_class, nan=0.5)

    binary_predictions = (predictions_positive_class >= 0.5).astype(int)

    accuracy = accuracy_score(true_labels, binary_predictions)
    precision = precision_score(true_labels, binary_predictions)
    recall = recall_score(true_labels, binary_predictions)
    f1 = f1_score(true_labels, binary_predictions)

    print(f"{model_name} - Accuracy: {accuracy:.4f}")
    print(f"{model_name} - Precision: {precision:.4f}")
    print(f"{model_name} - Recall: {recall:.4f}")
    print(f"{model_name} - F1-score: {f1:.4f}")

    conf_matrix = confusion_matrix(true_labels, binary_predictions)
    conf_matrix_df = pd.DataFrame(conf_matrix, columns=['Predicted 0', 'Predicted 1'], index=['Actual 0', 'Actual 1'])
    print("\nConfusion Matrix:")
    print(conf_matrix_df)

    fpr, tpr, thresholds = roc_curve(true_labels, predictions_positive_class)
    auc = roc_auc_score(true_labels, predictions_positive_class)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, label=f'{model_name}, AUC = {auc:.4f}')
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.legend()
    plt.show()
