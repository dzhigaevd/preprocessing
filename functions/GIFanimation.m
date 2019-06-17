function GIFanimation(imFileName, figureHandle, timeDelay, noOfImages, index)
    currentFrame = getframe(figureHandle);
    imAni = frame2im(currentFrame);
    [A,map] = rgb2ind(imAni,256);
    if index == 1;
        imwrite(A,map,imFileName,'gif','LoopCount',Inf,'DelayTime',timeDelay);
    elseif index == noOfImages
        imwrite(A,map,imFileName,'gif','WriteMode','append','DelayTime',3);
    else
        imwrite(A,map,imFileName,'gif','WriteMode','append','DelayTime',timeDelay);
    end
end