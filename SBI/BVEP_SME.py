#!/usr/bin/env python3
"""
@author: meysamhashemi  INS Marseille

"""
import os
import sys
import numpy as np
import re
import matplotlib.pyplot as plt

######################################################
def true_node_binary(true_labels_nodes):
    true_nodes=np.zeros((len(true_labels_nodes)))
    for i in np.arange(np.r_[len(true_labels_nodes)]):
        if str(true_labels_nodes[i])==str('EZ'):
                      true_nodes[i]=2.
        elif  str(true_labels_nodes[i])==str('PZ'):
                      true_nodes[i]=1.
        else:
                       true_nodes[i]=0.  
    return  true_nodes       
######################################################
def estimated_node_binary(est_labels_nodes):
    estimated_nodes=np.zeros((len(est_labels_nodes)))
    for i in np.arange(np.r_[len(est_labels_nodes)]):
        if str(est_labels_nodes[i])==str('EZ'):
                      estimated_nodes[i]=2.
        elif  str(est_labels_nodes[i])==str('PZ'):
                      estimated_nodes[i]=1.  
        else:
                       estimated_nodes[i]=0.  
    return  estimated_nodes       
######################################################
######################################################
def SME_labels(eta_true, eta_posterior, eta_c, delta_eta):
    true_labels_nodes=[]
    est_labels_nodes=[]
    eta_est=np.mean(eta_posterior, axis=1)
   
    if eta_est.shape[0]==eta_true.shape[0]:
        for i, idx in enumerate(eta_true):
                if  eta_c<=eta_true[i]<eta_c+delta_eta:
                    true_labels_nodes.append('EZ')
                elif  eta_c-delta_eta<=eta_true[i]<eta_c:
                    true_labels_nodes.append('PZ')
                else: 
                    true_labels_nodes.append('HZ')


                if  eta_c<=eta_est[i]<eta_c+delta_eta:
                    est_labels_nodes.append('EZ')
                elif  eta_c-delta_eta<=eta_est[i]<eta_c:
                    est_labels_nodes.append('PZ')
                else: 
                    est_labels_nodes.append('HZ')
            
    return true_labels_nodes, est_labels_nodes
######################################################

def SEM_violinplot(nn, eta_true, ez_idx, pz_idx, eta_posterior, eta_c, delta_eta ):
    parts= plt.violinplot(eta_posterior, widths=0.7, showmeans=True, showextrema=True);
    plt.plot(np.r_[0:nn]+1,eta_true ,'o', color='k', alpha=0.9, markersize=4)
    plt.axhline(y=eta_c, linewidth=.8, color = 'r', linestyle='--')
    plt.axhline(y=eta_c-delta_eta, linewidth=.8, color = 'y', linestyle='--')
    plt.yticks(fontsize=14) 
    plt.xticks(fontsize=14) 
    #plt.xticks(np.r_[1:nn+1], np.r_[1:nn+1], rotation=90, fontsize=14)  
    #plt.xticks(np.arange(1,nn+1, step=2),np.arange(1, nn+1, step=2), fontsize=12, rotation=0)
    plt.ylabel(' Posterior ' +r'${(\eta_i)}$', fontsize=22);  
    plt.xlabel('Brain nodes', fontsize=22); 

    for pc in parts['bodies'][0:nn]:
        pc.set_facecolor('g')
        pc.set_edgecolor('g')
        pc.set_alpha(0.5)
    i = 0
    while i < len(ez_idx):
        for pc in parts['bodies'][ez_idx[i]:ez_idx[i]+1]:
            pc.set_facecolor('r')
            pc.set_edgecolor('r')
            pc.set_alpha(0.8)
        i += 1

    j = 0
    while j < len(pz_idx):
        for pc in parts['bodies'][pz_idx[j]:pz_idx[j]+1]:
            pc.set_facecolor('y')
            pc.set_edgecolor('y')
            pc.set_alpha(0.8)
        j += 1
    plt.tight_layout()
    
######################################################
def plot_zscore_shrinkage(nodes, eta_true, eta_est_mu, eta_est_std, prior_std):
    z_score_eta=z_score(eta_true, eta_est_mu, eta_est_std)
    colors= np.random.rand(z_score_eta.shape[0])
    plt.scatter(shrinkage([prior_std]*nodes.shape[0], eta_est_std), z_score_eta ,s=120, c='blue')
    plt.xlabel("Posterior shrinkages", fontsize=24)
    plt.ylabel("Posterior z-scores", fontsize=24)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.axis((-.1,1.1,-2,20))  
######################################################
def plot_confusion_matrix(cm, target_names,  cmap=None, normalize=False):
    import numpy as np
    import itertools
    import matplotlib.pyplot as plt

    accuracy = np.trace(cm) / float(np.sum(cm))
    misclass = 1 - accuracy

    if cmap is None:
        cmap = plt.get_cmap('Blues')

    #plt.figure(figsize=(8, 6))
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    #plt.title(title)
    plt.colorbar()

    if target_names is not None:
        tick_marks = np.arange(len(target_names))
        plt.xticks(tick_marks, target_names, rotation=45)
        plt.yticks(tick_marks, target_names, rotation=45)

    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]


    thresh = cm.max() / 1.5 if normalize else cm.max() / 2
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        if normalize:
            plt.text(j, i, "{:0.4f}".format(cm[i, j]),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")
        else:
            plt.text(j, i, "{:,}".format(cm[i, j]),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black", fontsize=24)

    #plt.tight_layout()
    plt.title('Confusion matrix', fontsize=28)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.ylabel('True nodes', fontsize=28)
    plt.xlabel('Predicted nodes\naccuracy={:0.4f}; misclass={:0.4f}'.format(accuracy, misclass), fontsize=18)
    plt.xlabel('Predicted nodes', fontsize=28)


######################################################
def plot_roc_curve(fpr, tpr):
    import matplotlib.pyplot as plt
    plt.style.use('seaborn')
    lw=2
    roc_auc=metrics.auc(fpr, tpr)
    plt.plot(fpr, tpr, color='orange',  lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='darkblue', linestyle='--')
    plt.xlabel('False Positive Rate', fontsize=18)
    plt.ylabel('True Positive Rate', fontsize=18)
    plt.title('ROC for amortized inference of SEEG', fontsize=14)
    plt.legend(loc="lower right")
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.label=('ROC curve (area = %0.2f)' % roc_auc)
    plt.legend(loc="lower right")
    plt.show()
######################################################   

from sklearn import metrics
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
from sklearn.preprocessing import label_binarize
from itertools import cycle

######################################################

def SME_accuracy(true_labels_nodes, est_labels_nodes):
    node_classes=['HZ', 'PZ', 'EZ']
    n_classes =3  
    
    true_nodes=true_node_binary(true_labels_nodes)
    estimated_nodes=estimated_node_binary(est_labels_nodes)
  
    y_true=label_binarize(true_nodes, classes=[0, 1, 2])
    y_pred=label_binarize(estimated_nodes, classes=[0, 1, 2])

    
    fpr = dict()
    tpr = dict()
    roc_auc = dict()


    for i in range(n_classes):
        fpr[i], tpr[i], _ = metrics.roc_curve(y_true[:, i], y_pred[:, i])
        roc_auc[i] = metrics.auc(fpr[i], tpr[i])

    fpr["micro"], tpr["micro"], _ = metrics.roc_curve(y_true.ravel(), y_pred.ravel())
    roc_auc["micro"] = metrics.auc(fpr["micro"], tpr["micro"])    
    
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(n_classes):
        mean_tpr += interp(all_fpr, fpr[i], tpr[i])

    # Finally average it and compute AUC
    mean_tpr /= n_classes

    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = metrics.auc(fpr["macro"], tpr["macro"])
    
    return roc_auc["micro"], roc_auc["macro"]

######################################################

def _binary_clf_curve(y_true, y_score):
    """
    Calculate true and false positives per binary classification
    threshold (can be used for roc curve or precision/recall curve);
    the calcuation makes the assumption that the positive case
    will always be labeled as 1

    Parameters
    ----------
    y_true : 1d ndarray, shape = [n_samples]
        True targets/labels of binary classification

    y_score : 1d ndarray, shape = [n_samples]
        Estimated probabilities or scores

    Returns
    -------
    tps : 1d ndarray
        True positives counts, index i records the number
        of positive samples that got assigned a
        score >= thresholds[i].
        The total number of positive samples is equal to
        tps[-1] (thus false negatives are given by tps[-1] - tps)

    fps : 1d ndarray
        False positives counts, index i records the number
        of negative samples that got assigned a
        score >= thresholds[i].
        The total number of negative samples is equal to
        fps[-1] (thus true negatives are given by fps[-1] - fps)

    thresholds : 1d ndarray
        Predicted score sorted in decreasing order

    References
    ----------
    Github: scikit-learn _binary_clf_curve
    - https://github.com/scikit-learn/scikit-learn/blob/ab93d65/sklearn/metrics/ranking.py#L263
    """

    # sort predicted scores in descending order
    # and also reorder corresponding truth values
    desc_score_indices = np.argsort(y_score)[::-1]
    y_score = y_score[desc_score_indices]
    y_true = y_true[desc_score_indices]

    # y_score typically consists of tied values. Here we extract
    # the indices associated with the distinct values. We also
    # concatenate a value for the end of the curve
    distinct_indices = np.where(np.diff(y_score))[0]
    end = np.array([y_true.size - 1])
    threshold_indices = np.hstack((distinct_indices, end))

    thresholds = y_score[threshold_indices]
    tps = np.cumsum(y_true)[threshold_indices]

    # (1 + threshold_indices) = the number of positives
    # at each index, thus number of data points minus true
    # positives = false positives
    fps = (1 + threshold_indices) - tps
    return tps, fps, thresholds

######################################################

def _roc_auc_score(y_true, y_score):
    """
    Compute Area Under the Curve (AUC) from prediction scores

    Parameters
    ----------
    y_true : 1d ndarray, shape = [n_samples]
        True targets/labels of binary classification

    y_score : 1d ndarray, shape = [n_samples]
        Estimated probabilities or scores

    Returns
    -------
    auc : float
    """

    # ensure the target is binary
    if np.unique(y_true).size != 2:
        raise ValueError('Only two class should be present in y_true. ROC AUC score '
                         'is not defined in that case.')
    
    tps, fps, _ = _binary_clf_curve(y_true, y_score)

    # convert count to rate
    tpr = tps / tps[-1]
    fpr = fps / fps[-1]

    # compute AUC using the trapezoidal rule;
    # appending an extra 0 is just to ensure the length matches
    zero = np.array([0])
    tpr_diff = np.hstack((np.diff(tpr), zero))
    fpr_diff = np.hstack((np.diff(fpr), zero))
    auc = np.dot(tpr, fpr_diff) + np.dot(tpr_diff, fpr_diff) / 2
    return auc    
######################################################    