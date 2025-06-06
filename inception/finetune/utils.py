""" helper function

author baiyu
"""
import os
import sys
import re
import datetime

import numpy

import torch
from torch.optim.lr_scheduler import _LRScheduler
import torchvision
import torchvision.transforms as transforms
from torch.utils.data import DataLoader
from torch.utils.data import Dataset

import json
from jsonpath import jsonpath

from PIL import Image

import glob


def get_network(args):
    """ return given network
    """

    if args.net == 'vgg16':
        from models.vgg import vgg16_bn
        net = vgg16_bn()
    elif args.net == 'vgg13':
        from models.vgg import vgg13_bn
        net = vgg13_bn()
    elif args.net == 'vgg11':
        from models.vgg import vgg11_bn
        net = vgg11_bn()
    elif args.net == 'vgg19':
        from models.vgg import vgg19_bn
        net = vgg19_bn()
    elif args.net == 'densenet121':
        from models.densenet import densenet121
        net = densenet121()
    elif args.net == 'densenet161':
        from models.densenet import densenet161
        net = densenet161()
    elif args.net == 'densenet169':
        from models.densenet import densenet169
        net = densenet169()
    elif args.net == 'densenet201':
        from models.densenet import densenet201
        net = densenet201()
    elif args.net == 'googlenet':
        from models.googlenet import googlenet
        net = googlenet()
    elif args.net == 'inceptionv3':
        from models.inceptionv3 import inceptionv3
        net = inceptionv3()
    elif args.net == 'inceptionv4':
        from inceptionv4 import inceptionv4
        net = inceptionv4()
    elif args.net == 'inceptionresnetv2':
        from models.inceptionv4 import inception_resnet_v2
        net = inception_resnet_v2()
    elif args.net == 'xception':
        from models.xception import xception
        net = xception()
    elif args.net == 'resnet18':
        from models.resnet import resnet18
        net = resnet18()
    elif args.net == 'resnet34':
        from models.resnet import resnet34
        net = resnet34()
    elif args.net == 'resnet50':
        from models.resnet import resnet50
        net = resnet50()
    elif args.net == 'resnet101':
        from models.resnet import resnet101
        net = resnet101()
    elif args.net == 'resnet152':
        from models.resnet import resnet152
        net = resnet152()
    elif args.net == 'preactresnet18':
        from models.preactresnet import preactresnet18
        net = preactresnet18()
    elif args.net == 'preactresnet34':
        from models.preactresnet import preactresnet34
        net = preactresnet34()
    elif args.net == 'preactresnet50':
        from models.preactresnet import preactresnet50
        net = preactresnet50()
    elif args.net == 'preactresnet101':
        from models.preactresnet import preactresnet101
        net = preactresnet101()
    elif args.net == 'preactresnet152':
        from models.preactresnet import preactresnet152
        net = preactresnet152()
    elif args.net == 'resnext50':
        from models.resnext import resnext50
        net = resnext50()
    elif args.net == 'resnext101':
        from models.resnext import resnext101
        net = resnext101()
    elif args.net == 'resnext152':
        from models.resnext import resnext152
        net = resnext152()
    elif args.net == 'shufflenet':
        from models.shufflenet import shufflenet
        net = shufflenet()
    elif args.net == 'shufflenetv2':
        from models.shufflenetv2 import shufflenetv2
        net = shufflenetv2()
    elif args.net == 'squeezenet':
        from models.squeezenet import squeezenet
        net = squeezenet()
    elif args.net == 'mobilenet':
        from models.mobilenet import mobilenet
        net = mobilenet()
    elif args.net == 'mobilenetv2':
        from models.mobilenetv2 import mobilenetv2
        net = mobilenetv2()
    elif args.net == 'nasnet':
        from models.nasnet import nasnet
        net = nasnet()
    elif args.net == 'attention56':
        from models.attention import attention56
        net = attention56()
    elif args.net == 'attention92':
        from models.attention import attention92
        net = attention92()
    elif args.net == 'seresnet18':
        from models.senet import seresnet18
        net = seresnet18()
    elif args.net == 'seresnet34':
        from models.senet import seresnet34
        net = seresnet34()
    elif args.net == 'seresnet50':
        from models.senet import seresnet50
        net = seresnet50()
    elif args.net == 'seresnet101':
        from models.senet import seresnet101
        net = seresnet101()
    elif args.net == 'seresnet152':
        from models.senet import seresnet152
        net = seresnet152()
    elif args.net == 'wideresnet':
        from models.wideresidual import wideresnet
        net = wideresnet()
    elif args.net == 'stochasticdepth18':
        from models.stochasticdepth import stochastic_depth_resnet18
        net = stochastic_depth_resnet18()
    elif args.net == 'stochasticdepth34':
        from models.stochasticdepth import stochastic_depth_resnet34
        net = stochastic_depth_resnet34()
    elif args.net == 'stochasticdepth50':
        from models.stochasticdepth import stochastic_depth_resnet50
        net = stochastic_depth_resnet50()
    elif args.net == 'stochasticdepth101':
        from models.stochasticdepth import stochastic_depth_resnet101
        net = stochastic_depth_resnet101()

    else:
        print('the network name you have entered is not supported yet')
        sys.exit()

    if args.gpu: #use_gpu
        net = net.cuda()

    return net

class Mydatasetpro(Dataset):
# 类初始化
    def __init__(self, img_paths, labels, transform):
        self.imgs = img_paths
        # self.ims = img
        self.labels = labels
        # self.las = las
        self.transforms = transform
# 进行切片
    def __getitem__(self, index):                #根据给出的索引进行切片，并对其进行数据处理转换成Tensor，返回成Tensor
        img = self.imgs[index]
        # ims = self.ims[index]
        label = self.labels[index]
        # lal = self.las[index]
        #pil_img = Image.open(img)                 #pip install pillow
        pil_img = numpy.load(img)
        np_img = numpy.nan_to_num(pil_img)
        data = self.transforms(np_img)
        return data, label
# 返回长度
    def __len__(self):
        return len(self.imgs)

class Mydatasettest(Dataset):
# 类初始化
    def __init__(self, img_paths, labels, ims, transform):
        self.imgs = img_paths
        # self.ims = img
        self.labels = labels
        self.ims = ims
        # self.las = las
        self.transforms = transform
# 进行切片
    def __getitem__(self, index):                #根据给出的索引进行切片，并对其进行数据处理转换成Tensor，返回成Tensor
        img = self.imgs[index]
        # ims = self.ims[index]
        label = self.labels[index]
        ims = self.ims[index]
        # lal = self.las[index]
        #pil_img = Image.open(img)                 #pip install pillow
        pil_img = numpy.load(img)
        np_img = numpy.nan_to_num(pil_img)
        data = self.transforms(np_img)
        return data, label, ims
# 返回长度
    def __len__(self):
        return len(self.imgs)

class Mydatasetpre(Dataset):
# 类初始化
    def __init__(self, img_paths, ims, transform):
        self.imgs = img_paths
        self.ims = ims
        #self.labels = labels
        # self.las = las
        self.transforms = transform
# 进行切片
    def __getitem__(self, index):                #根据给出的索引进行切片，并对其进行数据处理转换成Tensor，返回成Tensor
        img = self.imgs[index]
        ims = self.ims[index]
        #label = self.labels[index]
        # lal = self.las[index]
        #pil_img = Image.open(img)                 #pip install pillow
        pil_img = numpy.load(img)
        np_img = numpy.nan_to_num(pil_img)
        data = self.transforms(np_img)
        return data, ims
# 返回长度
    def __len__(self):
        return len(self.imgs)

def get_training_dataloader(mean, std, batch_size=16, num_workers=2, shuffle=True):
    """ return training dataloader
    Args:
        mean: mean of cifar100 training dataset
        std: std of cifar100 training dataset
        path: path to cifar100 training python dataset
        batch_size: dataloader batchsize
        num_workers: dataloader num_works
        shuffle: whether to shuffle
    Returns: train_data_loader:torch dataloader object
    """

    labelpath = "inception/finetune/data/train_l"
    labelfiles = os.listdir(labelpath)
    # labelpath30 = "/root/autodl-tmp/whb/NA12878/label_chr20_d30"
    # labelfiles30= os.listdir(labelpath30)
    imgpath = "inception/finetune/data/train_i"
    imgfiles= os.listdir(imgpath)
    # imgpath30 = "/root/autodl-tmp/whb/NA12878/images_chr20_d30"
    # imgfiles30= os.listdir(imgpath30)

    imgfiles = list(filter(file_filter, imgfiles))
    labelfiles = list(filter(j_filter, labelfiles))
    # imgfiles30 = list(filter(file_filter, imgfiles30))

    all_imgs_path = []
    for imgf in imgfiles:
        all_imgs_path.append(imgpath+"/"+imgf)
    # for imgf in imgfiles30:
    #     all_imgs_path.append(imgpath30+"/"+imgf)

    # print(all_imgs_path)

    all_labels = []
    for lf in labelfiles:
        with open(labelpath+"/"+lf) as f:
            data = json.load(f)
            # print(data)
            label = jsonpath(data, "$..snp")
            # print(label)
            if label[0] == True:
                all_labels.append(0)
            else:
                all_labels.append(1)
    # for lf in labelfiles30:
    #     with open(labelpath30+"/"+lf) as f:
    #         data = json.load(f)
    #         # print(data)
    #         label = jsonpath(data, "$..snp")
    #         # print(label)
    #         if label[0] == True:
    #             all_labels.append(0)
    #         else:
    #             all_labels.append(1)

    # print(all_labels)

    index = numpy.random.permutation(len(all_imgs_path))

    train_imgs = numpy.array(all_imgs_path)[index]
    train_labels = numpy.array(all_labels)[index]


    transform_train = transforms.Compose([
        #transforms.ToPILImage(),
        # transforms.RandomCrop(32, padding=4),
        # transforms.RandomHorizontalFlip(),
        # transforms.RandomRotation(15),
        #transforms.Resize((32,32)),
        transforms.ToTensor()
        #transforms.Normalize(mean, std)
    ])
    #cifar100_training = CIFAR100Train(path, transform=transform_train)
    # cifar100_training = torchvision.datasets.CIFAR100(root='./data', train=True, download=True, transform=transform_train)
    cifar100_training = Mydatasetpro(train_imgs, train_labels, transform=transform_train) 
    cifar100_training_loader = DataLoader(
        cifar100_training, shuffle=shuffle, num_workers=num_workers, batch_size=batch_size)

    return cifar100_training_loader

def get_val_dataloader(mean, std, batch_size=16, num_workers=2, shuffle=True):
    """ return training dataloader
    Args:
        mean: mean of cifar100 test dataset
        std: std of cifar100 test dataset
        path: path to cifar100 test python dataset
        batch_size: dataloader batchsize
        num_workers: dataloader num_works
        shuffle: whether to shuffle
    Returns: cifar100_test_loader:torch dataloader object
    """
    labelpath = "inception/finetune/data/val_l"
    labelfiles= os.listdir(labelpath)
    # labelpath30 = "/root/autodl-tmp/pytorch-cifar100/data/val_l30"
    # labelfiles30= os.listdir(labelpath30)
    
    imgpath = "inception/finetune/data/val_i"
    imgfiles= os.listdir(imgpath)
    # imgpath30 = "/root/autodl-tmp/pytorch-cifar100/data/val_i30"
    # imgfiles30= os.listdir(imgpath30) 

    imgfiles = list(filter(file_filter, imgfiles))
    labelfiles = list(filter(j_filter, labelfiles))
    # imgfiles30 = list(filter(file_filter, imgfiles30))

    all_imgs_path = []
    for imgf in imgfiles:
        all_imgs_path.append(imgpath+"/"+imgf)
    # for imgf in imgfiles30:
    #     all_imgs_path.append(imgpath30+"/"+imgf)

    # print(all_imgs_path)

    all_labels = []
    for lf in labelfiles:
        with open(labelpath+"/"+lf) as f:
            data = json.load(f)
            # print(data)
            label = jsonpath(data, "$..snp")
            # print(label)
            if label[0] == True:
                all_labels.append(0)
            else:
                all_labels.append(1)
    # for lf in labelfiles30:
    #     with open(labelpath30+"/"+lf) as f:
    #         data = json.load(f)
    #         # print(data)
    #         label = jsonpath(data, "$..snp")
    #         # print(label)
    #         if label[0] == True:
    #             all_labels.append(0)
    #         else:
    #             all_labels.append(1)

    index = numpy.random.permutation(len(all_imgs_path))

    test_imgs = numpy.array(all_imgs_path)[index]
    test_labels = numpy.array(all_labels)[index]

    transform_test = transforms.Compose([
        
        #transforms.Resize((32,32)),
        transforms.ToTensor()
        #transforms.Normalize(mean, std)
    ])
    #cifar100_test = CIFAR100Test(path, transform=transform_test)
    # cifar100_test = torchvision.datasets.CIFAR100(root='./data', train=False, download=True, transform=transform_test)
    cifar100_test = Mydatasetpro(test_imgs, test_labels, transform=transform_test)
    cifar100_test_loader = DataLoader(
        cifar100_test, shuffle=shuffle, num_workers=num_workers, batch_size=batch_size)

    return cifar100_test_loader

def get_pred_dataloader(batch_size=16, num_workers=2, shuffle=False):
    """ return training dataloader
    Args:
        mean: mean of cifar100 test dataset
        std: std of cifar100 test dataset
        path: path to cifar100 test python dataset
        batch_size: dataloader batchsize
        num_workers: dataloader num_works
        shuffle: whether to shuffle
    Returns: cifar100_test_loader:torch dataloader object
    """
    #labelpath = "/root/autodl-tmp/SNPTools/label_d30"
    #labelfiles= os.listdir(labelpath)
    imgpath = "/root/autodl-tmp/SNPTools/Z.mays/images_5"
    imgfiles= os.listdir(imgpath)

    imgfiles = list(filter(file_filter, imgfiles))

    all_imgs_path = []
    for imgf in imgfiles:
        all_imgs_path.append(imgpath+"/"+imgf)

    # print(all_imgs_path)


    #index = numpy.random.permutation(len(all_imgs_path))

    #test_imgs = numpy.array(all_imgs_path)[index]
    #test_ims = numpy.array(imgfiles)[index]

    test_imgs = numpy.array(all_imgs_path)
    test_ims = numpy.array(imgfiles)
 
    transform_test = transforms.Compose([
        
        #transforms.Resize((32,32)),
        transforms.ToTensor()
        #transforms.Normalize(mean, std)
    ])
    #cifar100_test = CIFAR100Test(path, transform=transform_test)
    # cifar100_test = torchvision.datasets.CIFAR100(root='./data', train=False, download=True, transform=transform_test)
    cifar100_test = Mydatasetpre(test_imgs, test_ims, transform=transform_test)
    cifar100_test_loader = DataLoader(
        cifar100_test, shuffle=shuffle, num_workers=num_workers, batch_size=batch_size)

    return cifar100_test_loader


def file_filter(f):
    if f[-4:] in ['.npy']:
        return True
    else:
        return False

def j_filter(f):
    if f[-5:] in ['.json']:
        return True
    else:
        return False

def get_test_dataloader(mean, std, batch_size=16, num_workers=2, shuffle=True):
    """ return training dataloader
    Args:
        mean: mean of cifar100 test dataset
        std: std of cifar100 test dataset
        path: path to cifar100 test python dataset
        batch_size: dataloader batchsize
        num_workers: dataloader num_works
        shuffle: whether to shuffle
    Returns: cifar100_test_loader:torch dataloader object
    """
    labelpath = "/root/autodl-tmp/SNPTools/label_d30"
    labelfiles= os.listdir(labelpath)
    imgpath = "/root/autodl-tmp/SNPTools/images_d30"
    imgfiles= os.listdir(imgpath)

    imgfiles = list(filter(file_filter, imgfiles))

    all_imgs_path = []
    for imgf in imgfiles:
        all_imgs_path.append(imgpath+"/"+imgf)

    # print(all_imgs_path)

    all_labels = []
    for lf in labelfiles:
        with open(labelpath+"/"+lf) as f:
            data = json.load(f)
            # print(data)
            label = jsonpath(data, "$..snp")
            # print(label)
            if label[0] == True:
                all_labels.append(0)
            else:
                all_labels.append(1)

    # print(all_labels)

    index = numpy.random.permutation(len(all_imgs_path))

    test_imgs = numpy.array(all_imgs_path)[index]
    test_labels = numpy.array(all_labels)[index]
    test_ims = numpy.array(imgfiles)[index]

    transform_test = transforms.Compose([
        
        #transforms.Resize((32,32)),
        transforms.ToTensor()
        #transforms.Normalize(mean, std)
    ])
    #cifar100_test = CIFAR100Test(path, transform=transform_test)
    # cifar100_test = torchvision.datasets.CIFAR100(root='./data', train=False, download=True, transform=transform_test)
    cifar100_test = Mydatasettest(test_imgs, test_labels, test_ims, transform=transform_test)
    cifar100_test_loader = DataLoader(
        cifar100_test, shuffle=shuffle, num_workers=num_workers, batch_size=batch_size)

    return cifar100_test_loader
def compute_mean_std(cifar100_dataset):
    """compute the mean and std of cifar100 dataset
    Args:
        cifar100_training_dataset or cifar100_test_dataset
        witch derived from class torch.utils.data

    Returns:
        a tuple contains mean, std value of entire dataset
    """

    data_r = numpy.dstack([cifar100_dataset[i][1][:, :, 0] for i in range(len(cifar100_dataset))])
    data_g = numpy.dstack([cifar100_dataset[i][1][:, :, 1] for i in range(len(cifar100_dataset))])
    data_b = numpy.dstack([cifar100_dataset[i][1][:, :, 2] for i in range(len(cifar100_dataset))])
    mean = numpy.mean(data_r), numpy.mean(data_g), numpy.mean(data_b)
    std = numpy.std(data_r), numpy.std(data_g), numpy.std(data_b)

    return mean, std

class WarmUpLR(_LRScheduler):
    """warmup_training learning rate scheduler
    Args:
        optimizer: optimzier(e.g. SGD)
        total_iters: totoal_iters of warmup phase
    """
    def __init__(self, optimizer, total_iters, last_epoch=-1):

        self.total_iters = total_iters
        super().__init__(optimizer, last_epoch)

    def get_lr(self):
        """we will use the first m batches, and set the learning
        rate to base_lr * m / total_iters
        """
        return [base_lr * self.last_epoch / (self.total_iters + 1e-8) for base_lr in self.base_lrs]


def most_recent_folder(net_weights, fmt):
    """
        return most recent created folder under net_weights
        if no none-empty folder were found, return empty folder
    """
    # get subfolders in net_weights
    folders = os.listdir(net_weights)

    # filter out empty folders
    folders = [f for f in folders if len(os.listdir(os.path.join(net_weights, f)))]
    if len(folders) == 0:
        return ''

    # sort folders by folder created time
    folders = sorted(folders, key=lambda f: datetime.datetime.strptime(f, fmt))
    return folders[-1]

def most_recent_weights(weights_folder):
    """
        return most recent created weights file
        if folder is empty return empty string
    """
    weight_files = os.listdir(weights_folder)
    if len(weights_folder) == 0:
        return ''

    regex_str = r'([A-Za-z0-9]+)-([0-9]+)-(regular|best)'

    # sort files by epoch
    weight_files = sorted(weight_files, key=lambda w: int(re.search(regex_str, w).groups()[1]))

    return weight_files[-1]

def last_epoch(weights_folder):
    weight_file = most_recent_weights(weights_folder)
    if not weight_file:
       raise Exception('no recent weights were found')
    resume_epoch = int(weight_file.split('-')[1])

    return resume_epoch

def best_acc_weights(weights_folder):
    """
        return the best acc .pth file in given folder, if no
        best acc weights file were found, return empty string
    """
    files = os.listdir(weights_folder)
    if len(files) == 0:
        return ''

    regex_str = r'([A-Za-z0-9]+)-([0-9]+)-(regular|best)'
    best_files = [w for w in files if re.search(regex_str, w).groups()[2] == 'best']
    if len(best_files) == 0:
        return ''

    best_files = sorted(best_files, key=lambda w: int(re.search(regex_str, w).groups()[1]))
    return best_files[-1]
