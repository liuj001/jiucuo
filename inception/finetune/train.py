# train.py
#!/usr/bin/env	python3

""" train network using pytorch

author baiyu
"""

import os
import sys
import argparse
import time
from datetime import datetime

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import torchvision
import torchvision.transforms as transforms

from torch.utils.data import DataLoader
# from torch.utils.tensorboard import SummaryWriter

#from conf import settings
from utils import get_network, get_training_dataloader, get_val_dataloader, WarmUpLR, \
    most_recent_folder, most_recent_weights, last_epoch, best_acc_weights

def train(epoch):

    start = time.time()
    net.train()
    for batch_index, (images, labels) in enumerate(cifar100_training_loader):

        if args.gpu:
            #labels = labels.cuda()
            #labels = labels.tolist() 
            #if labels[0] > 1:
              #labels[0] = 1
            #elif labels[1] > 1:
              #labels[1] = 1
            #labels = torch.tensor(labels)
            labels = labels.cuda()
            images = images.cuda()

        optimizer.zero_grad()
        output = net(images)
        outputs = torch.nan_to_num(output)
        
        m=nn.Sigmoid()
        outputs_s=m(outputs)

        #print('outputs: {outputs}, labels: {labels}'.format(
        #  outputs=outputs_s,
        #    labels=labels
        #))

        loss = loss_function(outputs_s, labels)
        loss.backward()
        optimizer.step()

        n_iter = (epoch - 1) * len(cifar100_training_loader) + batch_index + 1

        last_layer = list(net.children())[-1]
        # for name, para in last_layer.named_parameters():
        #     if 'weight' in name:
        #         writer.add_scalar('LastLayerGradients/grad_norm2_weights', para.grad.norm(), n_iter)
        #     if 'bias' in name:
        #         writer.add_scalar('LastLayerGradients/grad_norm2_bias', para.grad.norm(), n_iter)

        print('Training Epoch: {epoch} [{trained_samples}/{total_samples}]\tLoss: {:0.4f}\tLR: {:0.8f}'.format(
            loss.item(),
            optimizer.param_groups[0]['lr'],
            epoch=epoch,
            trained_samples=batch_index * args.b + len(images),
            total_samples=len(cifar100_training_loader.dataset)
        ))

        #update training loss for each iteration
        #writer.add_scalar('Train/loss', loss.item(), n_iter)

        if epoch <= args.warm:
            warmup_scheduler.step()

    for name, param in net.named_parameters():
        layer, attr = os.path.splitext(name)
        attr = attr[1:]
        # writer.add_histogram("{}/{}".format(layer, attr), param, epoch)

    finish = time.time()

    print('epoch {} training time consumed: {:.2f}s'.format(epoch, finish - start))

@torch.no_grad()
def eval_training(epoch=0, tb=True):

    start = time.time()
    net.eval()

    test_loss = 0.0 # cost function error
    correct = 0.0

    for (images, labels) in cifar100_test_loader:

        if args.gpu:
            images = images.cuda()
            #labels = labels.tolist() 
            #if labels[0] > 1:
              #labels[0] = 1
            #elif labels[1] > 1:
              #labels[1] = 1
            #labels = torch.tensor(labels)
            labels = labels.cuda()

        output = net(images)
        outputs = torch.nan_to_num(output)

        m=nn.Sigmoid()
        outputs_s=m(outputs)

        #print('outputs: {outputs}, labels: {labels}'.format(
        #    outputs=outputs_s,
        #    labels=labels
        #))

        loss = loss_function(outputs_s, labels)

        test_loss += loss.item()
        _, preds = outputs_s.max(1)
        correct += preds.eq(labels).sum()

    finish = time.time()
    #if args.gpu:
        #print('GPU INFO.....')
        # print(torch.cuda.memory_summary(), end='')
    print('Evaluating Network.....')
    print('Test set: Epoch: {}, Average loss: {:.4f}, Accuracy: {:.4f}, Time consumed:{:.2f}s'.format(
        epoch,
        test_loss / len(cifar100_test_loader.dataset),
        correct.float() / len(cifar100_test_loader.dataset),
        finish - start
    ))
    print()

    #add informations to tensorboard
    # if tb:
        # writer.add_scalar('Test/Average loss', test_loss / len(cifar100_test_loader.dataset), epoch)
        # writer.add_scalar('Test/Accuracy', correct.float() / len(cifar100_test_loader.dataset), epoch)

    return correct.float() / len(cifar100_test_loader.dataset)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-net', type=str, default='inceptionv4', help='net type')
    parser.add_argument('-gpu', action='store_true', default=True, help='use gpu or not')
    parser.add_argument('-b', type=int, default=8, help='batch size for dataloader')
    parser.add_argument('-warm', type=int, default=2, help='warm up training phase')
    parser.add_argument('-lr', type=float, default=0.0001, help='initial learning rate')
    parser.add_argument('-resume', action='store_true', default=False, help='resume training')
    args = parser.parse_args()

    net = get_network(args)
    CIFAR100_TRAIN_MEAN = (0.5070751592371323, 0.48654887331495095, 0.4409178433670343)
    CIFAR100_TRAIN_STD = (0.2673342858792401, 0.2564384629170883, 0.27615047132568404)
    #directory to save weights file
    CHECKPOINT_PATH = 'inception/finetune/checkpoint'

    #total training epoches
    EPOCH = 200
    MILESTONES = [60, 120, 160]

    #EPOCH = 10
    #MILESTONES = [3, 7]

    #initial learning rate
    #INIT_LR = 0.1

    DATE_FORMAT = '%A_%d_%B_%Y_%Hh_%Mm_%Ss'
    #time of we run the script
    TIME_NOW = datetime.now().strftime(DATE_FORMAT)

    #tensorboard log dir
    LOG_DIR = 'runs'

    #save weights file per SAVE_EPOCH epoch
    SAVE_EPOCH = 10
    #SAVE_EPOCH = 1
    #data preprocessing:
    cifar100_training_loader = get_training_dataloader(
        CIFAR100_TRAIN_MEAN,
        CIFAR100_TRAIN_STD,
        num_workers=4,
        batch_size=args.b,
        shuffle=True
    )

    cifar100_test_loader = get_val_dataloader(
        CIFAR100_TRAIN_MEAN,
        CIFAR100_TRAIN_STD,
        num_workers=4,
        batch_size=args.b,
        shuffle=True
    )

    loss_function = nn.CrossEntropyLoss()
    #loss_function = nn.BCELoss()
    optimizer = optim.SGD(net.parameters(), lr=args.lr, momentum=0.9, weight_decay=5e-4)
    train_scheduler = optim.lr_scheduler.MultiStepLR(optimizer, milestones=MILESTONES, gamma=0.2) #learning rate decay
    iter_per_epoch = len(cifar100_training_loader)
    warmup_scheduler = WarmUpLR(optimizer, iter_per_epoch * args.warm)

    if args.resume:
        recent_folder = most_recent_folder(os.path.join(CHECKPOINT_PATH, args.net), fmt=settings.DATE_FORMAT)
        if not recent_folder:
            raise Exception('no recent folder were found')

        checkpoint_path = os.path.join(CHECKPOINT_PATH, args.net, recent_folder)

    else:
        checkpoint_path = os.path.join(CHECKPOINT_PATH, args.net, TIME_NOW)

    #use tensorboard
    if not os.path.exists(LOG_DIR):
        os.mkdir(LOG_DIR)

    #since tensorboard can't overwrite old values
    #so the only way is to create a new tensorboard log
    # writer = SummaryWriter(log_dir=os.path.join(
            # settings.LOG_DIR, args.net, settings.TIME_NOW))
    #input_tensor = torch.Tensor(1, 3, 32, 32)
    input_tensor = torch.Tensor(1, 6, 32, 32)
    if args.gpu:
        input_tensor = input_tensor.cuda()
    # writer.add_graph(net, input_tensor)

    #create checkpoint folder to save model
    if not os.path.exists(checkpoint_path):
        os.makedirs(checkpoint_path)
    checkpoint_path = os.path.join(checkpoint_path, '{net}-{epoch}-{type}.pth')

    best_acc = 0.0
    if args.resume:
        best_weights = best_acc_weights(os.path.join(CHECKPOINT_PATH, args.net, recent_folder))
        if best_weights:
            weights_path = os.path.join(CHECKPOINT_PATH, args.net, recent_folder, best_weights)
            print('found best acc weights file:{}'.format(weights_path))
            print('load best training file to test acc...')
            net.load_state_dict(torch.load(weights_path))
            best_acc = eval_training(tb=False)
            print('best acc is {:0.2f}'.format(best_acc))

        recent_weights_file = most_recent_weights(os.path.join(CHECKPOINT_PATH, args.net, recent_folder))
        if not recent_weights_file:
            raise Exception('no recent weights file were found')
        weights_path = os.path.join(CHECKPOINT_PATH, args.net, recent_folder, recent_weights_file)
        print('loading weights file {} to resume training.....'.format(weights_path))
        net.load_state_dict(torch.load(weights_path))

        resume_epoch = last_epoch(os.path.join(CHECKPOINT_PATH, args.net, recent_folder))


    for epoch in range(1, EPOCH + 1):
        if epoch > args.warm:
            # train_scheduler.step(epoch)
            train_scheduler.step()

        if args.resume:
            if epoch <= resume_epoch:
                continue

        train(epoch)
        acc = eval_training(epoch)

        #start to save best performance model after learning rate decay to 0.01
        if epoch > MILESTONES[1] and best_acc < acc:
            weights_path = checkpoint_path.format(net=args.net, epoch=epoch, type='best')
            print('saving weights file to {}'.format(weights_path))
            torch.save(net.state_dict(), weights_path)
            best_acc = acc
            continue

        if not epoch % SAVE_EPOCH:
            weights_path = checkpoint_path.format(net=args.net, epoch=epoch, type='regular')
            print('saving weights file to {}'.format(weights_path))
            torch.save(net.state_dict(), weights_path)

    # writer.close()
