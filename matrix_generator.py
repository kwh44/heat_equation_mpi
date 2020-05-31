import random

height = 50
width = 50

with open('matrix.txt', 'w') as matrix:
  for i in range(height):
    for j in range(width):
      if i == 0 or i == (height - 1) or j == 0:
        matrix.write(str(100) + ' ')
      elif j == (width - 1):
        matrix.write(str(100))
      else:
        matrix.write(str(0) + ' ')
    matrix.write('\n')
