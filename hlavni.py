# importing the eel library  
import eel

# initializing the application
eel.init("static_web_folder")

@eel.expose
def give_result(a, b):
    eel.print_status(18)()
    return (b)

# starting the application
eel.start("index.html")