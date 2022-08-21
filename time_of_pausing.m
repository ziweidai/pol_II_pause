function dt = time_of_pausing(v0,k,b,TES,gene_length)
total_time_without_pausing = (gene_length+1000)/v0;
total_time_with_pausing = integral(@(x) (((x-TES).^2+b)/v0)./((x-TES).^2+k*b),...
    0,gene_length+1000);
dt = total_time_with_pausing - total_time_without_pausing;
end