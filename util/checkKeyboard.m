function [newdown sec]=checkKeyboard()
  global downmap; 
 
  [keyIsDown sec keyCode] = KbCheckM;
  
  % find the newly pressed keys
  newdown = keyCode & xor(keyCode,downmap);
  if(~newdown)
    downmap = keyCode;
    return
  end

  % save the key presses
  newdown = find(newdown);
 
  downmap = keyCode;
return
    