function output = filtering(f,x)

% output at one time instant

% 1) まずベクトル化（行ベクトルへ）
f = f(:).';     % 1 × numel(f)
x = x(:).';     % 1 × numel(x)

if length(f) >= length(x)

    output= sum(f(1:1:length(x)).*x(end:-1:1));
else
    output = sum(f.*x(end:-1:end-length(f)+1));

end
end