import eel

eel.init('static_web_folder')

def print_num(n):
    print('Got this from Javascript:', n)



eel.start('index.html', size=(400, 300), block=False)

n = eel.js_random()()  # This immediately returns the value
print('Got this from Javascript:', n)