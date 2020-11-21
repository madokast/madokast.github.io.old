import torch

x = torch.tensor(1.)
a = torch.tensor(1.,requires_grad=True)
b = torch.tensor(2.,requires_grad=True)
c = torch.tensor(3.,requires_grad=True)

y = a**2 *b + b*x +c
print(y)
print(y.grad_fn)

print(a.grad,b.grad,c.grad)

grads = torch.autograd.grad(y,[a,b,c])


print(grads)