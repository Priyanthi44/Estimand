% 
for i=1:10
optimiseParameters(.5 ,3)

end
% Rf(1:10) = parallel.FevalFuture;
% for idx = 1:10
%     f(idx) = parfeval(@magic,1,idx);
% end
% magicResults = cell(1,10);
% for idx = 1:10
%     [completedIdx,value] = fetchNext(f);
%     magicResults{completedIdx} = value;
%     fprintf('Got result with index: %d.\n', completedIdx);
% end

%   f(1:3) = parallel.FevalFuture;
%         
%         f(1)= parfeval(@optimiseParameters, 1, .05 , beta);
%         
%         f(2) = parfeval(@runDNN, 1, 3, 1,1,7,false);
%         
%         f(3)=parfeval(@runDNN, 1,4, 1,1,7,true);
%         results = cell(1,3);
%       
%         
%         for idx = 1:3
%             [completedIdx,value] = fetchNext(f);
%             results{completedIdx} = value;
%             fprintf('Got result with index: %d. the value of %d\n', completedIdx, value);
%         end

% function testPersist()
%  persistent n
%     if isempty(n)
%         n = 0;
%     end
%     n = n+1;
% 
% end