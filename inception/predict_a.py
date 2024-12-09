#predict.py
#!/usr/bin/env python3

import os
os.environ["CUDA_VISIBLE_DEVICES"]="0"

import argparse

from matplotlib import pyplot as plt

import torch
import torchvision.transforms as transforms
from torch.utils.data import DataLoader

import numpy as np

import re
import pandas as pd

#from inception.conf import settings
from inception.utils_a import get_pred_dataloader
from inception.utils_a import get_network
#from utils_a import get_network, get_pred_dataloader

# import shutil

weights = "inception/weight/inceptionv4.pth"

def find_snp(snp_pic, snp_dir, chr):
    #print("in find_snp,"+snp_pic)
    cifar100_pred_loader = get_pred_dataloader(
        imgpath=snp_pic,
        batch_size=16,
        num_workers=4,
        shuffle=False
    )
    #print("data load done===============================================================")

    net = get_network()
    #net.cuda()
    #print("get_network+++++++++++++++++++++++++++++++++++++++")
    net.load_state_dict(torch.load(weights))
    #print("weights load done===============================================================")
    net.cuda().eval()
    #net.eval()
    #print("net load done===============================================================")

    pred_label_list = []
    true_ims_list = []
    #print("in pre===============================================================")
    with torch.no_grad():
        # for n_iter, (image, ims) in enumerate(cifar100_test_loader):
        #print("pre===============================================================")
        for n_iter, (image, ims) in enumerate(cifar100_pred_loader):
            #print("iteration: {}\ttotal {} iterations".format(n_iter + 1, len(cifar100_pred_loader)))

            image = image.cuda()
            #label = label.cuda()
                
                
            #print('GPU INFO.....')
            #print(torch.cuda.memory_summary(), end='')


            output = net(image)
            #_, pred = output.topk(5, 1, largest=True, sorted=True)
            _, pred = output.topk(1, largest=True, sorted=True)

            # label = label.view(label.size(0), -1).expand_as(pred)
            

            true_ims_list.append(ims)
            #true_label_list.append(label.cpu().detach().numpy())
            pred_label_list.append(pred.cpu().detach().numpy())
            # true_lal_list.append(lal)
    
    #print(snp_dir + '/'+ chr + '_snp.txt'+'++++++++++++++++++++++++++++++++++++++++++++++++++')
    y_ims = np.concatenate(true_ims_list)
    #y_true = np.concatenate(true_label_list)
    y_pred = np.concatenate(pred_label_list)
    # y_lal = np.concatenate(true_lal_list)
    #print("ims: ", y_ims)
    #print("pred: ", y_pred)

    snp_t = []
    snp_c = []
    for i in range(len(y_pred)):
        if y_pred[i]==0:
           #pos = re.sub('[snp.npy]', '', y_ims[i])
           #pos = int(pos)
           j = 0
           while y_ims[i][j] != '-':
               j += 1
           cname = y_ims[i][:j]
           #print("cname: ", cname)
           j += 1
           pos_o = y_ims[i][j:]
           #print("pos_o: ", pos_o)
           pos = re.sub('[snp.npy]', '', pos_o)
           pos = int(pos)
           snp_t.append(pos) 
           snp_c.append(cname)
    #print("pos: ", snp_t)
    
    out_file = open(snp_dir + '/'+ chr + '_snp.txt','wt')
    #print(snp_dir + '/'+ chr + '_snp.txt'+'++++++++++++++++++++++++++++++++++++++++++++++++++')
    #str = '\n'
    #out_file.write(str.join(snp_t))
    for i in range(len(snp_t)):
        out_file.write('{:}\t{:}\n'.format(snp_c[i], snp_t[i]))
    out_file.close()