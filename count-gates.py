import argparse
import re

def count_instr(d, name, v=1):
    if name in d:
        cnt = d[name]
    else:
        cnt = 0
    cnt += v
    d[name] = cnt

def parse_file(path, timer):
    instr = dict()
    in_timer = False
    with open(path, 'r') as fp:
        for line in fp.readlines():
            if not line.startswith('#'):
                instruction = re.split(r''',?\s''', line)
                name = instruction[0]
                if name == 'start' and timer is not None and int(instruction[1]) == timer:
                    if in_timer:
                        raise RuntimeError(f'Timer {timer} is already running')
                    in_timer = True
                elif name == 'stop' and timer is not None and int(instruction[1]) == timer:
                    if not in_timer:
                        raise RuntimeError(f'Timer {timer} is not running')
                    in_timer = False
                elif (timer is None or in_timer):# and not name.startswith('v'):
                    amount = 1
                    #if name == 'gmuls':
                    #    amount = int(instruction[1])//3
                        
                    count_instr(instr, name, amount)
                #if name.startswith('v'):
                #    raise RuntimeError(f'Vectorized instruction {name} not supported')
    return instr


parser = argparse.ArgumentParser()
parser.add_argument('-a', dest='file', type=str, required=True)
parser.add_argument('--block', dest='block', type=int, required=False)

args = parser.parse_args()
instr = parse_file(args.file, args.block)

for (name, cnt) in instr.items():
    print(f'{name}: {cnt}')